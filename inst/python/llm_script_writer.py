#!/usr/bin/env python3
"""
Given a dataset_id and output_dir, an LLM agent generates a bash script that appends
a markdown section to mm_dataset_<id>.md, self-corrects until it runs
cleanly, then prints the .md path.

Usage:
    python llm_script_writer.py <dataset_id> <output_dir>
"""
import asyncio
import os
import subprocess
import sys
from datetime import datetime

from agents import Agent, Runner, function_tool, set_default_openai_client, set_tracing_disabled
from openai import AsyncOpenAI

BASE_URL     = "http://fgcz-c-056:8081/v1"
MODEL        = "Qwen3.6-27B-FP8"
WORK_DIR     = "."  # overridden in __main__ to output_dir passed by caller
MAX_TURNS    = 15
EXEC_TIMEOUT = 30

set_tracing_disabled(True)
set_default_openai_client(AsyncOpenAI(base_url=BASE_URL, api_key="dummy", timeout=60.0, max_retries=1))

_succeeded: set[str] = set()


@function_tool
def write_file(path: str, content: str) -> str:
    """Write content to a file. Path must be under WORK_DIR."""
    if not os.path.abspath(path).startswith(WORK_DIR):
        raise ValueError(f"Path outside allowed directory: {path}")
    with open(path, "w") as f:
        f.write(content)
    print(f"  write_file: {path}")
    return f"Written: {path}"


@function_tool
def read_file(path: str) -> str:
    """Read a file's content. Reads are unrestricted so the agent can access analysis scripts."""
    if not os.path.exists(path):
        return f"File not found: {path}"
    with open(path, "r") as f:
        content = f.read()
    print(f"  read_file: {path}")
    return content


@function_tool
def run_script(path: str) -> dict:
    """Execute a bash script. Returns exit_code, stdout, stderr."""
    if not os.path.abspath(path).startswith(WORK_DIR):
        raise ValueError(f"Path outside allowed directory: {path}")
    if path in _succeeded:
        return {"exit_code": 0, "stdout": "Already succeeded.", "stderr": ""}
    try:
        r = subprocess.run(
            ["bash", path],
            capture_output=True,
            text=True,
            timeout=EXEC_TIMEOUT,
            cwd=WORK_DIR,
        )
        if r.returncode == 0 and not r.stderr.strip():
            _succeeded.add(path)
        print(f"  run_script: {path} -> exit_code={r.returncode}")
        return {"exit_code": r.returncode, "stdout": r.stdout, "stderr": r.stderr}
    except subprocess.TimeoutExpired:
        print(f"  run_script: {path} -> timeout")
        return {"exit_code": -1, "stdout": "", "stderr": f"Killed: script exceeded {EXEC_TIMEOUT}s"}


INSTRUCTIONS = """\
You write bash scripts and run them to produce output.
A script succeeds when exit_code is 0 AND stderr is empty.
If it fails, fix the script and retry. Once it succeeds, verify the output with read_file.
"""

agent = Agent(
    name="script_writer",
    model=MODEL,
    tools=[write_file, run_script, read_file],
    instructions=INSTRUCTIONS,
)


def _append_separator(dataset_id: int) -> str:
    """Append a run separator to the .md and return its path."""
    path = os.path.join(WORK_DIR, f"mm_dataset_{dataset_id}.md")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "a") as f:
        f.write(f"\n===== dataset_{dataset_id} | {timestamp} =====\n\n")
    return path


async def run(dataset_id: int) -> str:
    output_path = _append_separator(dataset_id)
    script_path = os.path.join(WORK_DIR, f"mm_dataset_{dataset_id}.sh")

    task = (
        f"Write a bash script to {script_path} that appends a Hello World markdown section "
        f"to {output_path} using >>. The script must exit 0 with empty stderr."
    )

    result = await Runner.run(
        starting_agent=agent,
        input=task,
        max_turns=MAX_TURNS,
    )

    if not os.path.exists(output_path):
        raise RuntimeError(f"Agent finished but output .md not found: {output_path}")

    return output_path


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python llm_script_writer.py <dataset_id> <output_dir>", file=sys.stderr)
        sys.exit(1)

    dataset_id = int(sys.argv[1])
    WORK_DIR = os.path.abspath(sys.argv[2])
    os.makedirs(WORK_DIR, exist_ok=True)
    path = asyncio.run(run(dataset_id))
    print(path)

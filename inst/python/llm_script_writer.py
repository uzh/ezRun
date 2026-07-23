#!/usr/bin/env python3
"""
LLM agent for Methods section generation.

Reads analysis scripts and logs, then appends a Methods paragraph per tool
directly to methods.txt.

Usage:
    python llm_script_writer.py --output <path> --identity-file <path> --task-file <path>
                                 [--scripts path ...] [--logs path ...]
"""
import argparse
import asyncio
import os
from datetime import datetime

from agents import Agent, Runner, function_tool, set_default_openai_client, set_tracing_disabled
from openai import AsyncOpenAI

BASE_URL  = "http://fgcz-c-056:8000/v1"
MODEL     = "DeepSeek-V4-Flash-DSpark"
OUTPUT_FILE = "methods.txt"  # overridden in __main__
WORK_DIR    = "."            # overridden in __main__ to dirname(OUTPUT_FILE)
MAX_TURNS = 25

set_tracing_disabled(True)
set_default_openai_client(AsyncOpenAI(base_url=BASE_URL, api_key="dummy", timeout=60.0, max_retries=1))


@function_tool
def append_file(path: str, content: str) -> str:
    """Append content to a file. Path must be under WORK_DIR."""
    if not os.path.abspath(path).startswith(WORK_DIR):
        raise ValueError(f"Path outside allowed directory: {path}")
    with open(path, "a") as f:
        f.write(content)
    print(f"  append_file: {path}")
    return f"Appended: {path}"


@function_tool
def read_file(path: str) -> str:
    """Read a file's content. Reads are unrestricted so the agent can access analysis scripts."""
    if not os.path.exists(path):
        return f"File not found: {path}"
    with open(path, "r") as f:
        content = f.read()
    print(f"  read_file: {path}")
    return content


def _output_path() -> str:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(OUTPUT_FILE, "a") as f:
        f.write(f"\n===== {timestamp} =====\n\n")
    return OUTPUT_FILE


async def run(script_paths: list[str], log_paths: list[str],
              identity: str, task_override: str) -> str:
    output_path = _output_path()

    effective_agent = Agent(
        name="methods_agent",
        model=MODEL,
        tools=[read_file, append_file],
        instructions=identity,
    )

    scripts_block = "\n".join(script_paths) if script_paths else "(none provided)"
    logs_block    = "\n".join(log_paths)    if log_paths    else "(none provided)"
    task_intro    = task_override

    task = f"""\
Read each of the following analysis scripts and log files using read_file.
Then append a Methods section to: {output_path}

{task_intro}

Use append_file to write the Methods text to {output_path}. Do not overwrite the file.

Analysis scripts:
{scripts_block}

Log files:
{logs_block}
"""

    size_before = os.path.getsize(output_path)

    result = await Runner.run(
        starting_agent=effective_agent,
        input=task,
        max_turns=MAX_TURNS,
    )

    if os.path.getsize(output_path) <= size_before:
        # LLM returned text without calling append_file — write final_output directly
        fallback = (result.final_output or "").strip()
        if fallback:
            with open(output_path, "a") as f:
                f.write(fallback + "\n")
            print(f"  fallback write: agent output captured from final_output")
        else:
            raise RuntimeError(f"Agent produced no output: {output_path} was not appended to")

    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",   required=True)
    parser.add_argument("--scripts",  nargs="*", default=[])
    parser.add_argument("--logs",     nargs="*", default=[])
    parser.add_argument("--identity-file", required=True)
    parser.add_argument("--task-file",     required=True)
    args = parser.parse_args()

    OUTPUT_FILE = os.path.abspath(args.output)
    WORK_DIR    = os.path.dirname(OUTPUT_FILE)
    os.makedirs(WORK_DIR, exist_ok=True)
    with open(args.identity_file) as f:
        identity = f.read()
    with open(args.task_file) as f:
        task = f.read()
    path = asyncio.run(run(args.scripts, args.logs, identity, task))
    print(path)

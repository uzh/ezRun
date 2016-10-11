/*!
 * this is probably no longer needed!!
 */


function popup(options) {
  var defaultOptions = {
    linkName: "",
    popup: true
  };
  
  if (typeof options.popup == 'undefined')
    options.popup = defaultOptions.popup;
  if (typeof options.linkName == 'undefined')
    alert('Missing link name.');
  
  var form = document.createElement('form');
  form.setAttribute('method', 'post');
  form.setAttribute('action', options.linkName);
  if (options.popup)
    form.setAttribute('target', '_blank');
  form.setAttribute('enctype', 'multipart/form-data');
  
  document.body.appendChild(form);
  form.submit();
  document.body.removeChild(form);
}

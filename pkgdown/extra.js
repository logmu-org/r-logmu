// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

// Re-render KaTeX inside Bootstrap popovers (used by pkgdown for footnotes)
// when they are shown, since KaTeX runs at page load before popover content
// is injected into the DOM.
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('[data-bs-toggle="popover"]').forEach(function (el) {
    el.addEventListener('shown.bs.popover', function () {
      var popoverId = el.getAttribute('aria-describedby');
      if (!popoverId) return;
      var popover = document.getElementById(popoverId);
      if (!popover) return;
      popover.querySelectorAll('.math.inline').forEach(function (span) {
        katex.render(span.textContent, span, { throwOnError: false });
      });
      popover.querySelectorAll('.math.display').forEach(function (span) {
        katex.render(span.textContent, span, { throwOnError: false, displayMode: true });
      });
    });
  });
});

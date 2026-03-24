// Break API reference h1 headings after each period
document.addEventListener('DOMContentLoaded', function() {
  // Only target h1 elements that contain "pterasoftware."
  const headings = document.querySelectorAll('h1');

  headings.forEach(function(h1) {
    const text = h1.textContent;

    // Only modify headings that look like module paths
    if (text.includes('pterasoftware.')) {
      // Add class for styling
      h1.classList.add('api-module-heading');

      // Get the headerlink anchor if it exists
      const headerlink = h1.querySelector('a.headerlink');

      // Split by period and rejoin with period + line break
      const parts = text.replace(/\s*¶\s*$/, '').split('.');

      // Clear the heading
      h1.innerHTML = '';

      // Add each part as a separate span with a line break after (except last)
      parts.forEach(function(part, index) {
        const span = document.createElement('span');
        span.textContent = part;
        h1.appendChild(span);

        if (index < parts.length - 1) {
          h1.appendChild(document.createTextNode('.'));
          h1.appendChild(document.createElement('br'));
        }
      });

      // Re-add the headerlink if it existed
      if (headerlink) {
        h1.appendChild(headerlink);
      }
    }
  });

  // Simplify right sidebar TOC entries by removing class name prefixes
  // e.g., "ClassName.method()" becomes "method()"
  const tocLinks = document.querySelectorAll('.toc-tree a');

  tocLinks.forEach(function(link) {
    const text = link.textContent;

    // Match patterns like "ClassName.method_name()" or "ClassName.attribute"
    // where ClassName starts with uppercase
    const match = text.match(/^[A-Z][A-Za-z0-9_]*\.(.+)$/);

    if (match) {
      link.textContent = match[1];
    }
  });
});

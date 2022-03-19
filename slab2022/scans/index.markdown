---
layout: home
---

<h1> Platinum layer target. March 2022</h1>

<a href="..">Up: slab target overview</a>

<h2>All scans analyzed to date </h2>
<ul>
{%- assign default_paths = site.pages | map: "path" -%}
{%- assign page_paths = site.header_pages | default: default_paths -%}
{%- for path in page_paths -%}
  {%- assign my_page = site.pages | where: "path", path | first -%}
  {%- if my_page.title -%}
  <li><a class="page-link" href="{{ my_page.url | relative_url }}">{{ my_page.title | escape }}</a></li>

  {%- endif -%}
{%- endfor -%}
</ul>

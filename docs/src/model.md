```@meta
EditURL = "https://github.com/sbfnk/bdbv-linelist-analysis/blob/main/MODEL.md"
```

```@eval
using Markdown
# Rewrite the GitHub-relative `[LIMITATIONS.md](LIMITATIONS.md)` link
# into a Vitepress-relative one so the docs build doesn't fail on a
# dead link while keeping the top-level MODEL.md readable on GitHub.
text = read(joinpath(@__DIR__, "..", "..", "MODEL.md"), String)
text = replace(text, "[LIMITATIONS.md](LIMITATIONS.md)" => "[Limitations](limitations)")
Markdown.parse(text)
```

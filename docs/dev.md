---
layout: default
title: Dev Notes
nav_order: 8
has_children: true
---

# Dev Notes

Internal implementation notes for contributors and maintainers. These pages are not required for using GeneCircuitry — see the [User Guide](../user-guide/) instead.

| Page                                                    | Contents                                                                |
| ------------------------------------------------------- | ----------------------------------------------------------------------- |
| [Controller Implementation](controller-implementation/) | Architecture diagrams, parallel processing flow, verification checklist |
| [Logging System](logging-system/)                       | `pipeline.log` and `error.log` formats, `log_step` / `log_error` API    |
| [Parallel Processing Fix](parallel-processing-fix/)     | Serialization fix for multiprocessing, pickle compatibility notes       |

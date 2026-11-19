# Discern explorer

Discern is a workbench for reviewing previously recorded evidence. Its main task is to open a saved isolate and connect its call to phenotype, provenance, and carrier evidence.

The interface uses a compact top navigation bar and a centered work surface. A low-chroma mineral ground, deep spruce controls, and separate bark-red, amber, and moss call colors distinguish actions from recorded findings. The label always carries the call meaning. System fonts keep the workspace available offline, and monospace is reserved for genomic keys and isolate IDs.

At laptop widths, the first viewport places the upload beside a preview of a real saved result. On the result page, the catalogue call and lab finding stay together, the saved score has its own panel, and the variant evidence is visible just below. Study details open on demand. Tables preserve header associations and scroll within their own region on narrow screens.

Navigation has a visible current-page state and keyboard focus. Local JavaScript adds drag-and-drop, file validation, and a submission state to the upload. The labelled native file picker stays available, and users can paste CSV text. Jinja and local assets render the entire interface offline.

[W3C table guidance](https://www.w3.org/WAI/tutorials/tables/) · [W3C focus visibility guidance](https://www.w3.org/WAI/WCAG22/Understanding/focus-not-obscured-minimum.html)

Component references include the [USWDS file input](https://designsystem.digital.gov/components/file-input/), [GOV.UK field errors](https://design-system.service.gov.uk/components/error-message/), [NHS summary list](https://service-manual.nhs.uk/design-system/components/summary-list), and [shadcn/ui cards](https://ui.shadcn.com/docs/components/base/card).

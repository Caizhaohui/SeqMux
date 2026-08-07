# NOTICE

SeqMux is an independent FASTQ demultiplexer.

Trimming ideas (e.g. BWA-style quality trim) follow well-known NGS practice also
used in tools such as Ultraplex / Cutadapt. SeqMux uses its **own sample barcode
table format** (`SampleNumber`, `Barcode1`, `Barcode2`) and does **not** implement
Ultraplex CSV compatibility.

References:

- Ultraplex: https://github.com/ulelab/ultraplex
- Cutadapt: https://github.com/marcelm/cutadapt

Both are distributed under the MIT License. SeqMux does not copy their source.

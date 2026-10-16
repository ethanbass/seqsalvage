#' The Sequencher project format
#'
#' Sequencher project files (`.SPF`) are undocumented. This page describes the
#' parts of the format that [read_spf()] relies on, as worked out from a
#' project saved by Sequencher 5.4 on a Mac. For the 150 reads in that project
#' whose original `.ab1` files survive, [read_spf()] reproduces the processed
#' traces, peak positions and quality scores exactly. The 14 projects published
#' with two studies (doi:10.1371/journal.pgen.1006176 and
#' [Zenodo record 5557681](https://zenodo.org/records/5557681)) have no original
#' files to compare against, but all 2,214 of their reads have one peak
#' position per base, all within the trace, and the called base is the
#' strongest trace channel at its peak for 94.8% of the 2.0 million A, C, G and
#' T calls. Projects from other versions may differ.
#'
#' All integers are big-endian and all offsets are absolute byte positions in
#' the file, counted from 0.
#'
#' # Overall layout
#'
#' The file is a heap of tagged records and data blocks, with no table of
#' contents at a fixed position. Records point at their data blocks by file
#' offset, so the blocks for one read are scattered through the file and do not
#' sit next to its record. Text such as read names, file paths, instrument names
#' and annotation labels is stored null-terminated. In projects saved on a Mac,
#' characters outside ASCII are Mac Roman: `c6` is ∆, as in `est2∆`. The
#' encoding of projects saved on Windows is not known; [read_spf()] and
#' [read_spf_contigs()] decode all text as Mac Roman.
#'
#' Records and blocks start with a 6-byte tag: a class byte, a version byte of
#' `30`, `31` or `32`, and `00 00 00 01`. In object tables and member lists the
#' class byte has bit `02` set when the version byte is `30` (`92`, `d2`) and
#' clear otherwise (`90`, `d0`). Both forms can occur in one file.
#'
#' # Read records
#'
#' Each imported trace has one record that starts with the 10 bytes
#' `c0 31 00 00 00 01 63 43 72 6d` (the last four spell `cCrm`), or with `32`
#' in place of `31`. The layout that follows is the same in both. Records
#' that start with `c2` instead appear to be reads deleted from the project:
#' their block tables point at space that has since been reused, and they are
#' skipped. Deleted reads can also keep the `c0` tag; see "Deleted records"
#' below.
#' [read_spf()] finds reads by searching for these tags. Relative to the start of
#' the tag:
#'
#' | Offset | Size | Contents |
#' |---|---|---|
#' | 44 | 4 | offset of the read's name, the file name it was imported under |
#' | 48 | 4 | length of the name, including the terminating null |
#' | 259 | 64 | block table: 8 pairs of (offset, length), 4 bytes each |
#'
#' The block table lists, in order:
#'
#' | Slot | Contents | Length field |
#' |---|---|---|
#' | 1 | base calls | number of bases |
#' | 2 | quality scores; offset and length are 0 if the read has none | number of bases |
#' | 3 | unused; offset and length are 0 | |
#' | 4 | peak positions, a `DnFx` block | compressed size |
#' | 5–8 | A, C, G and T trace channels, `DnFx` blocks | compressed size |
#'
#' The name shown in Sequencher can differ from the original file name: reads
#' can be renamed after import, and the record keeps only the current name.
#'
#' # Base calls and quality scores
#'
#' Base calls take one byte each, holding a 4-bit mask with A = 1, C = 2, G = 4
#' and T = 8. Ambiguity codes are the sums of these (for example R = A + G = 5),
#' and N is 15. [read_spf()] maps the masks to IUPAC letters and maps 0 to `-`.
#' In the published projects, the two bases of a two-base code are the two
#' strongest trace channels at that peak for 105 of 117 calls. These are
#' the calls as edited in Sequencher, so they can differ from the original file
#' where bases were inserted, deleted or changed.
#'
#' Quality scores take one byte each and are the original Phred scores.
#'
#' # `DnFx` blocks
#'
#' Peak positions and trace channels are stored as compressed series of 16-bit
#' integers. Each block starts with a 24-byte header:
#'
#' | Offset | Size | Contents |
#' |---|---|---|
#' | 0 | 4 | the ASCII tag `DnFx` |
#' | 4 | 4 | uncompressed size in bytes, two per value |
#' | 8 | 4 | size of the payload that follows the header |
#' | 12 | 4 | unknown; varies between blocks, possibly a checksum |
#' | 16 | 4 | `48 b7 ff bf` for peak positions, `68 b7 ff bf` for traces |
#' | 20 | 1 | 0 |
#' | 21 | 1 | order k of the difference predictor, 2 for peaks and 3–5 for traces |
#' | 22 | 2 | 0 |
#'
#' The payload starts with a chunk index:
#'
#' * bytes 0–1: the number of complete 256-value chunks, `floor(n / 256)`;
#' * bytes 2–3: `01 00` in every block seen;
#' * one 14-byte entry per chunk, holding the last five values of the chunk as
#'   16-bit integers, a flag byte (`0x80` or `0x00`, meaning unknown) and a
#'   3-byte position in the residual stream counted in nibbles.
#'
#' The index allows random access into a long trace. [read_spf()] does not
#' need it, because the residual stream that follows the index decodes from
#' start to finish without breaks at chunk boundaries.
#'
#' # Predictor
#'
#' Each value is predicted from the k values before it by a k-th order
#' difference: for k = 2 the prediction is `2 * x[i - 1] - x[i - 2]`, and in
#' general it is the sum over j = 1, ..., k of
#' `(-1)^(j + 1) * choose(k, j) * x[i - j]`. Values before the start of the
#' series count as 0. The stream stores the residual, the value minus its
#' prediction, so the series is recovered by summing the residuals cumulatively
#' k times.
#'
#' # Residual coding
#'
#' The residual stream is read four bits at a time, high nibble first.
#'
#' * A nibble from `0` to `7` or from `9` to `f` is a residual from -7 to 7,
#'   read as a signed 4-bit integer.
#' * The nibble `8` switches to bytes for one or more residuals. In each byte
#'   the top bit is set on the last residual of the run, after which the stream
#'   returns to nibbles. The low seven bits are a signed integer s, which stands
#'   for a residual of s + 7 if s > 0 and s - 7 if s < 0.
#' * If s = 0, the next four nibbles hold a signed 16-bit integer S, which
#'   stands for S + 7 or S - 7 in the same way. Large residuals at the start of
#'   a series are stored this way.
#'
#' For example, `3 d 0 7` is the residuals 3, -3, 0 and 7; `8 81` is 8 and
#' `8 ff` is -8; `8 7e 06 fb` is the run -9, 13, -12.
#'
#' # Object IDs
#'
#' Records refer to one another by 32-bit object ID rather than by offset.
#' Object tables are tagged `92` or `90`, then hold a 16-bit count and that
#' many pairs of (file offset, ID). Pages of member IDs (below) share the `90`
#' tag; [read_spf_contigs()] takes a block as an object table only if every
#' offset in it starts a record. Some IDs that records refer to are missing
#' from the tables as counted, but their (offset, ID) pairs still occur in the
#' file; [read_spf_contigs()] finds them by searching for the ID directly after
#' an offset that starts a record.
#'
#' # Contigs
#'
#' Contig records start with `c0 31 00 00 00 01 00 00 00 00`, or with `32` in
#' place of `31`. Relative to the start of the record:
#'
#' | Offset | Size | Contents |
#' |---|---|---|
#' | 28 | 1 | record type: 2 for a contig, 4 for a placement record |
#' | 48 | 1 | number of members |
#' | 49 | 4 | offset of the member list |
#' | 53 | 4 | unknown; 11 in every contig seen |
#' | 71 | 8 | offset and length of the consensus block |
#' | 205 | 8 | offset and length of the contig's name, null-terminated |
#'
#' A member list tagged `d2` or `d0` holds a 16-bit count and one 4-byte object
#' ID per member. A list tagged `52` or `50` is paged: its count is the number
#' of pages, followed by the file offset of each page, and each page is a block
#' tagged `92` or `90` with a count and up to 32 member IDs. Records of type 0
#' use the same layout; in some projects one of them, named after the project
#' file, lists the contigs.
#'
#' Each member is a placement record, which holds
#' the object ID of its read or reference sequence at offset 79. Reference
#' sequences (records tagged `cBAS`, two bases per byte) have no name of their
#' own; the placement record gives it at offset 205.
#'
#' The consensus block holds one byte per alignment column. The low four bits
#' are the base mask used for reads, and 0 where the column has no consensus
#' base; the meaning of the high four bits is unknown. It is the consensus as
#' last saved, not one recomputed from the reads. In the Sequencher 5.4 project,
#' every contig held only reads of the accession it was named after, plus the
#' reference. In the 14 published projects, every read that is not deleted
#' belongs to a contig, and each contig's member count at offset 48 matches its
#' member list. A contig record that starts with `c2`, like the deleted reads,
#' had no live reads and is skipped.
#'
#' # Deleted records
#'
#' Sequencher can leave the records of deleted reads and contigs in the file,
#' with their `c0` tags and data intact. [read_spf()] and [read_spf_contigs()]
#' count a record as deleted when no object table lists it, and return such
#' records only with `deleted = TRUE`. In 14 of the 16 projects examined, an
#' object table lists every read. One published project lists all but one, and
#' the unlisted read has the same traces as a listed one, so it is an old copy.
#' The last project lists 4 of its 30 reads. The other 26 decode cleanly; 24
#' carry names from a 2011 sequencing plate, and two have names overwritten by
#' later data. [read_spf()] names such reads `unnamed_read_<n>`.
#'
#' # Not yet decoded
#'
#' The position, orientation and gapped sequence of each read within its
#' contig, which the placement records also hold; annotation records (tagged
#' `c2 31 00 00`) such as the end trimming applied by Sequencher; and import
#' records tagged `<ID:`.
#'
#' @seealso [read_spf()], [spf_to_abif()], [write_abif()]
#' @name spf_format
NULL

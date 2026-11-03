# GDELT `QuadClass` against `EventRootCode`

Checked 2026-11-03 against raw GDELT 1.0 event files downloaded the same
day from <http://data.gdeltproject.org/events/> with
`download_gdelt(dest_dir = <temporary directory>)` (MD5-verified), read
with no row filters:

```r
read_gdelt(f, columns = c("EventRootCode", "EventCode", "QuadClass"),
           actor_types = NULL, cross_border = FALSE)
```

The files were deleted afterwards.

| File | Type | Rows |
|---|---|---:|
| `20200101.export.CSV.zip` | daily | 89,215 |
| `200601.zip` | monthly backfile | 314,884 |
| `1995.zip` | yearly backfile | 1,780,524 |

In all three files `EventRootCode` equals the first two characters of
`EventCode` in every row, and no row has a missing `QuadClass` or
`EventRootCode`.

## Cross-tabulations (`EventRootCode` x `QuadClass`, row counts)

Every root occurs with exactly one `QuadClass` in every file; empty cells
are zero.

| Root | 20200101 | 200601 | 1995 | `QuadClass` |
|---|---:|---:|---:|---:|
| 01 | 11,593 | 46,644 | 218,075 | 1 |
| 02 | 5,952 | 21,254 | 103,753 | 1 |
| 03 | 5,646 | 24,544 | 158,774 | 1 |
| 04 | 19,789 | 82,337 | 532,116 | 1 |
| 05 | 7,477 | 21,606 | 134,701 | 1 |
| 06 | 2,016 | 5,918 | 29,337 | 2 |
| 07 | 2,812 | 7,765 | 35,471 | 2 |
| 08 | 3,770 | 11,421 | 62,364 | 2 |
| 09 | 1,782 | 5,502 | 24,422 | 2 |
| 10 | 1,064 | 5,326 | 25,588 | 3 |
| 11 | 6,409 | 20,039 | 103,964 | 3 |
| 12 | 2,116 | 8,700 | 50,500 | 3 |
| 13 | 1,311 | 5,270 | 25,816 | 3 |
| 14 | 1,519 | 2,875 | 16,377 | 3 |
| 15 | 616 | 1,410 | 9,730 | 4 |
| 16 | 939 | 3,839 | 22,829 | 4 |
| 17 | 4,996 | 14,086 | 75,953 | 4 |
| 18 | 1,719 | 6,196 | 25,580 | 4 |
| 19 | 7,666 | 20,079 | 123,791 | 4 |
| 20 | 23 | 73 | 1,383 | 4 |

## Result

GDELT's `QuadClass` is a function of the root code, the same in all three
files:

| `QuadClass` | Label | Roots |
|---:|---|---|
| 1 | Verbal cooperation | 01–05 |
| 2 | Material cooperation | 06–09 |
| 3 | Verbal conflict | 10–14 |
| 4 | Material conflict | 15–20 |

This is the mapping `assign_quad()` already uses, so `cameo_lookup`'s
`QuadClass` matches GDELT's without change. The often-quoted alternative
(01–05 / 06–08 / 09–13 / 14–20) is not what GDELT 1.0 distributes.

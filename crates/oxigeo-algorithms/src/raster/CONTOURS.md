# Contour lines, from first principles

This explains what `generate_contours` (in `contour.rs`) computes and how,
starting from what a contour line is. Everything here describes the code as
it is; the worked example below is checked by the test
`test_contours_md_worked_example`.

## 1. What a contour line is

A digital elevation model (DEM), or any raster of a continuous quantity, is a
surface `z = f(x, y)` that we only know at regularly spaced sample points. A
**contour line** (or isoline) at level `L` is the set of points where the
surface is exactly `L`: `f(x, y) = L`. On a topographic map, the 100 m
contour is the curve you would walk along without going up or down, at a
height of 100 m.

Two facts about a continuous surface shape everything that follows:

- A contour line separates ground **above** `L` from ground **below** `L`.
  Walk from a point above to a point below, and you must cross it.
- A contour line never just stops. It either closes on itself (a ring
  around a hill or a hollow) or runs off the edge of the area we know about.

## 2. The grid

`generate_contours(data, width, height, config)` takes `width × height`
values in row-major order: the value at column `col`, row `row` is
`data[row * width + col]`. We treat that value as the surface's height at
the point `(x, y) = (col, row)`. So `x` grows to the right, `y` grows
**down** (row 0 is the top row), and grid points are 1 unit apart. Output
coordinates are in these units; a caller turns them into map coordinates.
(The georaster web app treats each value as the centre of a raster cell, so
it maps `(x, y)` to `(west + (x + 0.5) × cell width, north − (y + 0.5) ×
cell height)`.)

Four neighbouring grid points form a **cell**, the unit square with corners

```
(col, row)     tl ───── tr   (col + 1, row)
                │        │
                │  cell  │
                │        │
(col, row + 1) bl ───── br   (col + 1, row + 1)
```

A `width × height` grid has `(width − 1) × (height − 1)` cells. Each cell
edge is shared with at most one other cell: the cell on the other side of it,
unless the edge is on the border of the grid.

## 3. Which levels

`ContourConfig { interval, base, nodata }` asks for every level of the form
`base + k × interval` (for whole numbers `k`) that lies **strictly between**
the lowest and highest valid value in the grid. A level equal to the very
lowest or highest value is left out: it would only touch the surface at
isolated points, not cross it. Levels come out in ascending order, and so
does the output.

## 4. Above or below

For one level `L`, mark every grid point as **above** if its value is
`≥ L` and **below** if it is `< L`. (Points exactly at `L` count as above.
Something has to decide them, and this choice puts the line exactly
through such a point, see section 9.)

Because a contour separates above from below, it must cross every cell edge
whose two ends are marked differently, and no edge whose ends are marked the
same. (Strictly, a surface could dip below `L` and back between two points
that are both above. Between samples we have no information, so we assume
the surface is as simple as possible there.)

## 5. Where it crosses an edge: linear interpolation

Along an edge from value `v0` to value `v1`, assume the surface changes
linearly. It reaches `L` at the fraction

```
t = (L − v0) / (v1 − v0)
```

of the way along. On the top edge of the cell at `(col, row)` that is the
point `(col + t, row)`; on the left edge, `(col, row + t)`; and so on. If
`|v1 − v0|` is below `1e−15`, the edge is flat to within rounding and both
ends are a hair from `L`, so the formula is unreliable; the midpoint
`t = 0.5` is used instead.

The two cells that share an edge both compute its crossing from the same two
values, in the same direction, so they get exactly the same point. Section 7
relies on that.

## 6. Marching squares: the 16 cases

Each of a cell's four corners is above or below, so a cell has `2⁴ = 16`
possible patterns. The code numbers them with one bit per corner:
`tl = 8, tr = 4, bl = 2, br = 1`, a bit set when that corner is above. The
pattern fixes which edges are crossed, and the crossings are joined in pairs
by straight **segments**:

| Case | Corners above | Segments (edge to edge) |
|---|---|---|
| 0 | none | none |
| 1 | br | bottom – right |
| 2 | bl | left – bottom |
| 3 | bl, br | left – right |
| 4 | tr | top – right |
| 5 | tr, br | top – bottom |
| 6 | tr, bl | saddle, see below |
| 7 | tr, bl, br | top – left |
| 8 | tl | top – left |
| 9 | tl, br | saddle, see below |
| 10 | tl, bl | top – bottom |
| 11 | tl, bl, br | top – right |
| 12 | tl, tr | left – right |
| 13 | tl, tr, br | left – bottom |
| 14 | tl, tr, bl | bottom – right |
| 15 | all | none |

A case and its opposite (`k` and `15 − k`) cross the same edges, because
swapping above and below doesn't move the boundary between them. Every
crossed edge belongs to exactly one segment of the cell.

Running this over every cell is **marching squares**: march across the grid
one square at a time and emit each square's segments independently.

### Saddles (cases 6 and 9)

When diagonally opposite corners are above (tr and bl, or tl and br), all
four edges are crossed and there are two ways to pair the crossings:

Take case 9, tl and br above. Either the high ground runs diagonally from
tl to br through the middle of the cell, and the two lines cut off the low
corners; or tl and br are two separate high points, and the two lines cut
off the high corners:

| Case | Centre above (`≥ L`) | Centre below (`< L`) |
|---|---|---|
| 9: tl, br above | top – right, left – bottom: cuts off tr and bl; tl and br joined | top – left, bottom – right: cuts off tl and br; tl and br separate |
| 6: tr, bl above | top – left, bottom – right: cuts off tl and br; tr and bl joined | top – right, left – bottom: cuts off tr and bl; tr and bl separate |

The corners alone can't tell these apart; the surface inside the cell
decides. The code estimates the height at the cell's centre as the mean of
its four corners. That is exactly the centre value of the *bilinear*
surface, the smoothest surface through the four corners. If the centre is
`≥ L`, it is above, so the two above corners are joined through the middle;
otherwise the two below corners are.

This "centre decider" is the common practical choice. The exact answer for
a bilinear surface, the *asymptotic decider* (Nielson & Hamann, 1991),
compares `L` with the value at the surface's saddle point,
`(tl·br − tr·bl) / (tl + br − tr − bl)`, rather than at the centre. The two
can disagree when the level passes very close to the saddle. Either choice
gives valid, non-crossing lines; they differ only in which way a near-tie
connects.

## 7. Joining segments into lines

Marching squares gives a soup of short segments, each joining two points, in
no particular order: one or two per crossed cell. Users want whole lines, so the segments for
each level are joined end to end.

The key observation: **every segment end lies on a cell edge, and each edge
is shared by at most two cells.** So at most two segment ends meet at any
edge: one from the cell on each side. That makes joining a walk, with no
searching:

1. Give every edge an integer key. The horizontal edge from grid point
   `(col, row)` to `(col + 1, row)` is `2 × (row × width + col)`, and the
   vertical edge from `(col, row)` to `(col, row + 1)` is that plus 1. Each
   segment records the keys of the two edges its ends lie on.
2. Build an index from each edge key to the (one or two) segments ending
   there.
3. **Open lines.** A segment end whose edge has no second segment is the
   end of a line: it is on the grid border, or next to a cell skipped for
   nodata. Start at such an end and walk: leave the segment through its other
   edge, look up the segment on the far side of that edge, and continue until
   reaching an edge with nothing on the far side. That is the other end of
   the line.
4. **Closed lines.** Every segment not visited after that belongs to a ring
   (a line can only be open by ending somewhere). Start at any unvisited
   segment and walk until you are back at it. A closed line's last point
   repeats its first.

Each segment is visited once and each lookup is a hash-map access, so
joining takes time proportional to the number of segments.

Matching edges is also *exact*: two segments meet only if they share an
edge, which section 5 guarantees gives the same point. No floating-point
tolerance is involved.

### Why this replaced matching points

The previous joiner compared each new segment's end points (within `1e−9`)
with the ends of every line built so far, inserting at the front of a line
when it attached there, and then merged lines pairwise until none shared an
end. Each step is cheap, but on a large grid there are many segments and the
lines are long, so the total grew roughly with the square of the grid's
size: 36.6 s for the 2048 × 2048 benchmark in section 10, against 0.39 s
now. The two joiners produce the same lines; the test
`test_chaining_by_edge_matches_chaining_by_point` checks that on random grids
with holes.

## 8. Nodata

A value is treated as missing if it is NaN, infinite, or within `1e−10` of
`config.nodata`. Missing values are ignored when finding the range of
levels, and any cell with a missing corner is skipped entirely: we know
nothing about the surface there. Lines therefore end where they run into a
hole, as they do at the grid border.

## 9. Values exactly at a level

Integer DEMs with integer contour intervals often have grid values exactly
equal to a level. Such a point counts as above (section 4), so on an edge
from it to a below neighbour the crossing is at `t = 0`: the grid point
itself. Several edges meeting at that point can then all have their crossing
there, and two separate lines may **touch** at it. The joiner still follows
edges, not points, so each line continues through the cell across the edge
it arrived by, and lines never merge or swap at a touching point. The test
`test_chaining_uses_every_segment_once_when_values_equal_the_level` covers
this.

## 10. Putting it together

For each cell, once:

1. skip it if any corner is missing;
2. find the levels between its lowest and highest corner. A cell can only
   cross those. `levels` is sorted, so this is two binary searches;
3. for each such level, classify the corners, look up the case, interpolate
   the crossings, and append the segments to that level's list.

Then join each level's segments (section 7), and return a `ContourLine
{ level, points, is_closed }` for every line, in ascending order of level.
Within a level, open lines come first.

Visiting each cell once, rather than once per level, matters when there are
many levels: most cells cross none or one of them. Together with the edge
walk, contouring a synthetic terrain grid at 31 levels, natively and
single-threaded:

| | before | after |
|---|---|---|
| 256 × 256 | 0.50 s | 0.03 s |
| 1024 × 1024 | 8.6 s | 0.17 s |
| 2048 × 2048 | 36.6 s | 0.39 s |

These timings include writing the lines out as GeoJSON in the georaster web
app.

## 11. A worked example

Take a 3 × 3 grid with a peak in the middle and level `L = 3.5`
(`ContourConfig::new(3.5)` with the default base 0: 3.5 is the only multiple
of 3.5 strictly between the grid's minimum 1 and maximum 6):

```
         col 0   col 1   col 2
row 0      1       2       3
row 1      2       6       4
row 2      1       3       2
```

Above (`≥ 3.5`) are the 6 and the 4; everything else is below. There are
four cells:

- **Cell (0, 0)**, corners 1, 2 / 2, 6: only br is above, case 1,
  bottom – right. Bottom edge, 2 → 6: `t = (3.5 − 2) / 4 = 0.375`, point
  `(0.375, 1)`. Right edge, 2 → 6: `t = 0.375`, point `(1, 0.375)`.
- **Cell (1, 0)**, corners 2, 3 / 6, 4: bl and br above, case 3,
  left – right. Left edge, 2 → 6: `(1, 0.375)`, the same point as the cell
  beside it (section 5). Right edge, 3 → 4: `t = 0.5`, `(2, 0.5)`.
- **Cell (0, 1)**, corners 2, 6 / 1, 3: tr above, case 4, top – right. Top
  edge `(0.375, 1)`, shared with cell (0, 0). Right edge, 6 → 3:
  `t = (3.5 − 6) / (3 − 6) = 5/6`, point `(1, 1.833…)`.
- **Cell (1, 1)**, corners 6, 4 / 3, 2: tl and tr above, case 12,
  left – right. Left edge `(1, 1.833…)`, shared with cell (0, 1). Right
  edge, 4 → 2: `t = 0.25`, `(2, 1.25)`.

The right edges of cells (1, 0) and (1, 1) are on the grid border, so they
are the line's two ends. Walking from the first:

```
(2, 0.5) → (1, 0.375) → (0.375, 1) → (1, 1.833…) → (2, 1.25)
```

That's one open line of four segments, curving around the peak from the
right border back to the right border. That is all of the 3.5 contour this
grid can show: the 4 at `(2, 1)` is above too, so the ground stays above
3.5 out to the border between the two ends.

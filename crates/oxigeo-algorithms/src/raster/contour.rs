//! Contour line generation using the marching squares algorithm
//!
//! Extracts isolines (contour lines) from a 2D elevation grid at specified
//! intervals. The marching squares algorithm classifies each 2x2 cell by
//! comparing corner values to the contour level, then interpolates edge
//! crossings and chains them into polylines by walking from edge to edge.
//!
//! `CONTOURS.md`, next to this file, explains the algorithm from first
//! principles, with a worked example.

use crate::error::{AlgorithmError, Result};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// A point on a contour line
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContourPoint {
    pub x: f64,
    pub y: f64,
}

impl ContourPoint {
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// A single contour line at a specific elevation level
#[derive(Debug, Clone)]
pub struct ContourLine {
    pub level: f64,
    pub points: Vec<ContourPoint>,
    pub is_closed: bool,
}

/// Configuration for contour generation
#[derive(Debug, Clone)]
pub struct ContourConfig {
    /// Contour interval (must be > 0)
    pub interval: f64,
    /// Base level (default 0.0)
    pub base: f64,
    /// Nodata value to skip
    pub nodata: Option<f64>,
}

impl ContourConfig {
    /// Create a new contour configuration with the given interval.
    ///
    /// # Errors
    /// Returns an error if `interval` is not positive.
    pub fn new(interval: f64) -> Result<Self> {
        if interval <= 0.0 || !interval.is_finite() {
            return Err(AlgorithmError::InvalidParameter {
                parameter: "interval",
                message: format!("interval must be positive, got {interval}"),
            });
        }
        Ok(Self {
            interval,
            base: 0.0,
            nodata: None,
        })
    }

    #[must_use]
    pub fn with_base(mut self, base: f64) -> Self {
        self.base = base;
        self
    }

    #[must_use]
    pub fn with_nodata(mut self, nodata: f64) -> Self {
        self.nodata = Some(nodata);
        self
    }

    /// Compute the contour levels that span [min_val, max_val].
    fn compute_levels(&self, min_val: f64, max_val: f64) -> Vec<f64> {
        if min_val > max_val || !min_val.is_finite() || !max_val.is_finite() {
            return vec![];
        }

        let start = ((min_val - self.base) / self.interval).ceil();
        let end = ((max_val - self.base) / self.interval).floor();

        let start_i = start as i64;
        let end_i = end as i64;

        let mut levels = Vec::new();
        for i in start_i..=end_i {
            let level = self.base + (i as f64) * self.interval;
            if level > min_val && level < max_val {
                levels.push(level);
            }
        }
        levels
    }
}

/// Identifies the grid edge a contour point lies on: the edge between two
/// horizontally or vertically adjacent grid points.
///
/// Every point marching squares produces lies on a cell edge, and each edge
/// is shared by the (at most two) cells on either side of it, so two
/// segments that meet always share an edge, and computed their shared point
/// from the same two grid values. Matching segments by edge is exact,
/// unlike comparing their points' coordinates.
type EdgeKey = u64;

/// The edge from grid point `(col, row)` to `(col + 1, row)`.
#[inline]
fn horizontal_edge(col: usize, row: usize, width: usize) -> EdgeKey {
    ((row * width + col) as u64) << 1
}

/// The edge from grid point `(col, row)` to `(col, row + 1)`.
#[inline]
fn vertical_edge(col: usize, row: usize, width: usize) -> EdgeKey {
    (((row * width + col) as u64) << 1) | 1
}

/// A segment produced by marching squares for one cell, with the edge each
/// end lies on.
#[derive(Debug, Clone, Copy)]
struct Segment {
    p0: ContourPoint,
    p1: ContourPoint,
    e0: EdgeKey,
    e1: EdgeKey,
}

/// The segment from `a` to `b`, each a point and the edge it lies on.
#[inline]
fn segment(a: (EdgeKey, ContourPoint), b: (EdgeKey, ContourPoint)) -> Segment {
    Segment {
        p0: a.1,
        p1: b.1,
        e0: a.0,
        e1: b.0,
    }
}

/// Linearly interpolate between two values to find where `level` falls.
#[inline]
fn lerp_frac(v0: f64, v1: f64, level: f64) -> f64 {
    let denom = v1 - v0;
    if denom.abs() < 1e-15 {
        0.5
    } else {
        (level - v0) / denom
    }
}

/// Generate segments for a single cell using marching squares.
///
/// Corners ordered:
///   tl(0) --- tr(1)
///    |          |
///   bl(2) --- br(3)
///
/// Bit assignment: tl=bit3, tr=bit2, bl=bit1, br=bit0
#[allow(clippy::too_many_arguments)]
fn cell_segments(
    col: usize,
    row: usize,
    width: usize,
    tl: f64,
    tr: f64,
    bl: f64,
    br: f64,
    level: f64,
) -> Vec<Segment> {
    let case = ((tl >= level) as u8) << 3
        | ((tr >= level) as u8) << 2
        | ((bl >= level) as u8) << 1
        | (br >= level) as u8;

    // Where the level crosses each edge, by linear interpolation, and which
    // edge that is.
    // top edge: tl -> tr
    let top = || {
        let t = lerp_frac(tl, tr, level);
        (
            horizontal_edge(col, row, width),
            ContourPoint::new(col as f64 + t, row as f64),
        )
    };
    // bottom edge: bl -> br
    let bottom = || {
        let t = lerp_frac(bl, br, level);
        (
            horizontal_edge(col, row + 1, width),
            ContourPoint::new(col as f64 + t, row as f64 + 1.0),
        )
    };
    // left edge: tl -> bl
    let left = || {
        let t = lerp_frac(tl, bl, level);
        (
            vertical_edge(col, row, width),
            ContourPoint::new(col as f64, row as f64 + t),
        )
    };
    // right edge: tr -> br
    let right = || {
        let t = lerp_frac(tr, br, level);
        (
            vertical_edge(col + 1, row, width),
            ContourPoint::new(col as f64 + 1.0, row as f64 + t),
        )
    };

    match case {
        0 | 15 => vec![], // all below or all above
        1 => vec![segment(bottom(), right())],
        2 => vec![segment(left(), bottom())],
        3 => vec![segment(left(), right())],
        4 => vec![segment(top(), right())],
        5 => vec![segment(top(), bottom())],
        6 => {
            // Saddle point (TR and BL above): disambiguate using center value
            let center = (tl + tr + bl + br) * 0.25;
            if center >= level {
                // The centre joins TR to BL, cutting off TL and BR.
                vec![segment(top(), left()), segment(bottom(), right())]
            } else {
                // TR and BL are separate peaks.
                vec![segment(top(), right()), segment(left(), bottom())]
            }
        }
        7 => vec![segment(top(), left())],
        8 => vec![segment(top(), left())],
        9 => {
            // Saddle point (TL and BR above): disambiguate using center value
            let center = (tl + tr + bl + br) * 0.25;
            if center >= level {
                // The centre joins TL to BR, cutting off TR and BL.
                vec![segment(top(), right()), segment(left(), bottom())]
            } else {
                // TL and BR are separate peaks.
                vec![segment(top(), left()), segment(bottom(), right())]
            }
        }
        10 => vec![segment(top(), bottom())],
        11 => vec![segment(top(), right())],
        12 => vec![segment(left(), right())],
        13 => vec![segment(left(), bottom())],
        14 => vec![segment(bottom(), right())],
        _ => vec![], // unreachable for 4-bit
    }
}

/// Check if a value is considered nodata.
#[inline]
fn is_nodata(val: f64, nodata: Option<f64>) -> bool {
    if !val.is_finite() {
        return true;
    }
    if let Some(nd) = nodata {
        (val - nd).abs() < 1e-10
    } else {
        false
    }
}

/// Hashes an [`EdgeKey`] by Fibonacci multiplication. Edge keys are small,
/// distinct integers derived from grid positions, not attacker-chosen input,
/// so the standard library's DoS-resistant SipHash only costs time here.
#[derive(Default)]
struct EdgeHasher(u64);

impl Hasher for EdgeHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ u64::from(b)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }

    fn write_u64(&mut self, key: u64) {
        self.0 = key.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

type EdgeIndex = HashMap<EdgeKey, [usize; 2], BuildHasherDefault<EdgeHasher>>;

/// Marks an empty slot in the edge index.
const NO_SEGMENT: usize = usize::MAX;

/// Joins one level's segments into polylines.
///
/// A segment's two ends lie on two edges of its cell, and an edge's crossing
/// is shared with the segment of the cell on the other side of it, if there
/// is one (there isn't at the grid border, or next to a nodata cell). So each
/// edge carries at most two segment ends, and joining segments is a walk:
/// leave a segment through one of its edges, enter the segment on the other
/// side of that edge, and repeat until reaching an edge with no other
/// segment (an open line) or the starting segment (a closed one). Every
/// segment is visited once, so this is linear in the number of segments.
///
/// Open lines are walked first, from their ends, so each comes out whole;
/// every segment left after that belongs to a closed line. A closed line's
/// last point repeats its first.
fn chain_segments(segments: &[Segment]) -> Vec<(Vec<ContourPoint>, bool)> {
    // Which segments end on each edge.
    let mut at_edge = EdgeIndex::with_capacity_and_hasher(segments.len() * 2, Default::default());
    for (i, seg) in segments.iter().enumerate() {
        for edge in [seg.e0, seg.e1] {
            let ends = at_edge.entry(edge).or_insert([NO_SEGMENT; 2]);
            if ends[0] == NO_SEGMENT {
                ends[0] = i;
            } else if ends[1] == NO_SEGMENT {
                ends[1] = i;
            }
            // A third end on one edge is impossible for marching squares
            // output; if it ever happened that segment would simply end
            // there instead of joining.
        }
    }
    // The segment on the other side of `edge` from segment `from`, if any.
    let across = |edge: EdgeKey, from: usize| -> Option<usize> {
        let ends = at_edge.get(&edge)?;
        let next = if ends[0] == from { ends[1] } else { ends[0] };
        (next != NO_SEGMENT && next != from).then_some(next)
    };

    let mut visited = vec![false; segments.len()];
    let walk =
        |start: usize, start_edge: EdgeKey, visited: &mut [bool]| -> (Vec<ContourPoint>, bool) {
            let first = &segments[start];
            let mut points = vec![if first.e0 == start_edge {
                first.p0
            } else {
                first.p1
            }];
            let (mut current, mut entered) = (start, start_edge);
            loop {
                visited[current] = true;
                let seg = &segments[current];
                let (exit, point) = if seg.e0 == entered {
                    (seg.e1, seg.p1)
                } else {
                    (seg.e0, seg.p0)
                };
                points.push(point);
                match across(exit, current) {
                    Some(next) if next == start => return (points, true),
                    Some(next) if !visited[next] => {
                        current = next;
                        entered = exit;
                    }
                    _ => return (points, false),
                }
            }
        };

    let mut lines = Vec::new();
    // Open lines: start from each end that no other segment shares.
    for (i, seg) in segments.iter().enumerate() {
        if visited[i] {
            continue;
        }
        if across(seg.e0, i).is_none() {
            lines.push(walk(i, seg.e0, &mut visited));
        } else if across(seg.e1, i).is_none() {
            lines.push(walk(i, seg.e1, &mut visited));
        }
    }
    // Everything left is part of a closed line.
    for (i, seg) in segments.iter().enumerate() {
        if !visited[i] {
            lines.push(walk(i, seg.e0, &mut visited));
        }
    }
    lines
}

/// Generate contour lines from a 2D grid using marching squares.
///
/// # Arguments
/// * `data` - Row-major 2D grid of elevation values (width * height)
/// * `width` - Grid width
/// * `height` - Grid height
/// * `config` - Contour generation configuration
///
/// # Errors
/// Returns an error if `data.len() != width * height`.
///
/// # Returns
/// A `Vec<ContourLine>` sorted by level.
pub fn generate_contours(
    data: &[f64],
    width: usize,
    height: usize,
    config: &ContourConfig,
) -> Result<Vec<ContourLine>> {
    // Empty grids produce no contours
    if width == 0 || height == 0 {
        return Ok(vec![]);
    }

    if data.len() != width * height {
        return Err(AlgorithmError::InvalidDimensions {
            message: "data length must equal width * height",
            actual: data.len(),
            expected: width * height,
        });
    }

    // Need at least a 2x2 grid for marching squares
    if width < 2 || height < 2 {
        return Ok(vec![]);
    }

    // Find min/max, skipping nodata
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    for &v in data {
        if is_nodata(v, config.nodata) {
            continue;
        }
        if v < min_val {
            min_val = v;
        }
        if v > max_val {
            max_val = v;
        }
    }

    if !min_val.is_finite() || !max_val.is_finite() {
        return Ok(vec![]);
    }

    let levels = config.compute_levels(min_val, max_val);
    if levels.is_empty() {
        return Ok(vec![]);
    }

    // One pass over the cells. A cell only crosses the levels between its
    // lowest and highest corner (a level is crossed when some corner is at
    // or above it and some corner below it), and `levels` is ascending, so
    // each cell looks up just that range instead of every level re-reading
    // every cell. Each level's segments still come out in row-major cell
    // order.
    let mut segments_by_level: Vec<Vec<Segment>> = vec![Vec::new(); levels.len()];
    for row in 0..height - 1 {
        for col in 0..width - 1 {
            let tl = data[row * width + col];
            let tr = data[row * width + col + 1];
            let bl = data[(row + 1) * width + col];
            let br = data[(row + 1) * width + col + 1];

            // Skip cell if any corner is nodata
            if is_nodata(tl, config.nodata)
                || is_nodata(tr, config.nodata)
                || is_nodata(bl, config.nodata)
                || is_nodata(br, config.nodata)
            {
                continue;
            }

            let low = tl.min(tr).min(bl).min(br);
            let high = tl.max(tr).max(bl).max(br);
            let first = levels.partition_point(|&level| level <= low);
            let end = levels.partition_point(|&level| level <= high);
            for (index, &level) in levels.iter().enumerate().take(end).skip(first) {
                segments_by_level[index]
                    .extend_from_slice(&cell_segments(col, row, width, tl, tr, bl, br, level));
            }
        }
    }

    // Chain segments into polylines per level and collect results
    let mut contour_lines: Vec<ContourLine> = Vec::new();
    for (&level, segments) in levels.iter().zip(&segments_by_level) {
        for (points, is_closed) in chain_segments(segments) {
            if points.len() >= 2 {
                contour_lines.push(ContourLine {
                    level,
                    points,
                    is_closed,
                });
            }
        }
    }

    Ok(contour_lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contour_config_valid() {
        let config = ContourConfig::new(10.0);
        assert!(config.is_ok());
        let config = config.expect("should be ok");
        assert!((config.interval - 10.0).abs() < f64::EPSILON);
        assert!((config.base - 0.0).abs() < f64::EPSILON);
        assert!(config.nodata.is_none());
    }

    #[test]
    fn test_contour_config_zero_interval_error() {
        let result = ContourConfig::new(0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_contour_config_negative_interval_error() {
        let result = ContourConfig::new(-5.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_levels() {
        let config = ContourConfig::new(10.0).expect("valid config");
        let levels = config.compute_levels(5.0, 35.0);
        assert_eq!(levels, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_contour_flat_grid() {
        let data = vec![100.0; 16];
        let config = ContourConfig::new(10.0).expect("valid config");
        let contours = generate_contours(&data, 4, 4, &config).expect("should succeed");
        assert!(contours.is_empty());
    }

    #[test]
    fn test_contour_simple_slope() {
        // 4x4 grid with linear ramp from 0 to 30
        // Row 0: 0  10  20  30
        // Row 1: 0  10  20  30
        // Row 2: 0  10  20  30
        // Row 3: 0  10  20  30
        let data = vec![
            0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0,
            30.0,
        ];
        let config = ContourConfig::new(15.0).expect("valid config");
        let contours = generate_contours(&data, 4, 4, &config).expect("should succeed");

        // Should have contour at level 15.0 (between 0-30, strictly inside)
        assert!(!contours.is_empty());
        let levels: Vec<f64> = contours.iter().map(|c| c.level).collect();
        assert!(levels.contains(&15.0));

        // Each contour line should have at least 2 points
        for c in &contours {
            assert!(c.points.len() >= 2);
        }
    }

    #[test]
    fn test_contour_single_level() {
        // 3x3 grid:
        //  0   5  10
        //  0   5  10
        //  0   5  10
        let data = vec![0.0, 5.0, 10.0, 0.0, 5.0, 10.0, 0.0, 5.0, 10.0];
        let config = ContourConfig::new(3.0).expect("valid config");
        let contours = generate_contours(&data, 3, 3, &config).expect("should succeed");

        // Should produce contours at 3, 6, 9
        let levels: Vec<f64> = contours.iter().map(|c| c.level).collect();
        assert!(levels.contains(&3.0));
        assert!(levels.contains(&6.0));
        assert!(levels.contains(&9.0));
    }

    #[test]
    fn test_contour_nodata_handling() {
        // 3x3 grid with nodata in center
        let nodata = -9999.0;
        let data = vec![0.0, 5.0, 10.0, 0.0, nodata, 10.0, 0.0, 5.0, 10.0];
        let config = ContourConfig::new(4.0)
            .expect("valid config")
            .with_nodata(nodata);
        let contours = generate_contours(&data, 3, 3, &config).expect("should succeed");

        // Contours should still be generated, but cells touching nodata are skipped.
        // The center cell is nodata, so the 4 cells that touch it are all skipped.
        // Only corners that don't touch center could produce segments, but in a 3x3
        // grid all 4 cells touch the center — so no contours produced.
        // Actually, let's verify: cells are (0,0),(1,0),(0,1),(1,1).
        // Cell (0,0): corners are data[0],data[1],data[3],data[4]=nodata → skip
        // Cell (1,0): corners are data[1],data[2],data[4]=nodata,data[5] → skip
        // Cell (0,1): corners are data[3],data[4]=nodata,data[6],data[7] → skip
        // Cell (1,1): corners are data[4]=nodata,data[5],data[7],data[8] → skip
        // All cells skipped because they all touch center nodata.
        assert!(contours.is_empty());
    }

    #[test]
    fn test_contour_closed_loop() {
        // A "hill" grid: high value in center, low values on edges
        // 5x5 grid
        let data = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0, 0.0, 0.0, 5.0, 20.0, 5.0, 0.0, 0.0, 5.0,
            5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let config = ContourConfig::new(10.0).expect("valid config");
        let contours = generate_contours(&data, 5, 5, &config).expect("should succeed");

        // Should have a contour at level 10.0
        assert!(!contours.is_empty());
        let at_10: Vec<&ContourLine> = contours
            .iter()
            .filter(|c| (c.level - 10.0).abs() < 1e-10)
            .collect();
        assert!(!at_10.is_empty());

        // The contour around the central peak should be closed
        let has_closed = at_10.iter().any(|c| c.is_closed);
        assert!(has_closed, "expected a closed contour around the hill");
    }

    #[test]
    fn test_contour_saddle_point() {
        // Create a configuration that produces saddle points (cases 6/9).
        // Diagonal pattern: high in TL and BR, low in TR and BL (or vice versa).
        // 3x3 grid:
        //  10   0   10
        //   0  10    0
        //  10   0   10
        let data = vec![10.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0];
        let config = ContourConfig::new(4.0).expect("valid config");
        let contours = generate_contours(&data, 3, 3, &config).expect("should succeed");

        // Should produce contours (the saddle cells are disambiguated)
        assert!(!contours.is_empty());

        // All contour lines should have valid points
        for c in &contours {
            assert!(c.points.len() >= 2);
            for p in &c.points {
                assert!(p.x.is_finite() && p.y.is_finite());
            }
        }
    }

    /// Each line of a single-cell grid as its two endpoints, ordered, so
    /// segments can be compared regardless of direction.
    fn cell_lines(data: &[f64], level: f64) -> Vec<[(f64, f64); 2]> {
        let config = ContourConfig::new(level).expect("valid config");
        let contours = generate_contours(data, 2, 2, &config).expect("should succeed");
        let mut lines: Vec<[(f64, f64); 2]> = contours
            .iter()
            .filter(|c| c.level == level)
            .map(|c| {
                let a = c.points.first().expect("a point");
                let b = c.points.last().expect("a point");
                let (a, b) = ((a.x, a.y), (b.x, b.y));
                if a <= b { [a, b] } else { [b, a] }
            })
            .collect();
        lines.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
        lines
    }

    #[test]
    fn test_contour_east_west_slope_is_one_straight_line() {
        // The same ramp as `test_contour_simple_slope`, rising west to east:
        // the 15 contour is the straight line x = 1.5, from top to bottom.
        let data = vec![
            0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0,
            30.0,
        ];
        let config = ContourConfig::new(15.0).expect("valid config");
        let contours = generate_contours(&data, 4, 4, &config).expect("should succeed");
        let at_15: Vec<_> = contours.iter().filter(|c| c.level == 15.0).collect();
        assert_eq!(at_15.len(), 1, "one line: {at_15:?}");
        let line = at_15[0];
        assert!(
            line.points.iter().all(|p| (p.x - 1.5).abs() < 1e-12),
            "{line:?}"
        );
        let ys: Vec<f64> = line.points.iter().map(|p| p.y).collect();
        let (min_y, max_y) = ys
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), &y| (lo.min(y), hi.max(y)));
        assert_eq!((min_y, max_y), (0.0, 3.0));
    }

    #[test]
    fn test_contour_vertical_cells_cross_top_to_bottom() {
        // Cases 5 (right side above) and 10 (left side above).
        assert_eq!(
            cell_lines(&[0.0, 10.0, 0.0, 10.0], 5.0),
            vec![[(0.5, 0.0), (0.5, 1.0)]]
        );
        assert_eq!(
            cell_lines(&[10.0, 0.0, 10.0, 0.0], 5.0),
            vec![[(0.5, 0.0), (0.5, 1.0)]]
        );
    }

    #[test]
    fn test_contour_saddle_cells() {
        // Case 9: TL and BR above. With the centre (5) above the level, TL
        // and BR join and the line cuts off TR and BL; below it, TL and BR
        // are cut off instead.
        let tl_br = [10.0, 0.0, 0.0, 10.0];
        assert_eq!(
            cell_lines(&tl_br, 4.0),
            vec![[(0.0, 0.6), (0.4, 1.0)], [(0.6, 0.0), (1.0, 0.4)]]
        );
        assert_eq!(
            cell_lines(&tl_br, 6.0),
            vec![[(0.0, 0.4), (0.4, 0.0)], [(0.6, 1.0), (1.0, 0.6)]]
        );
        // Case 6: TR and BL above, the mirror image.
        let tr_bl = [0.0, 10.0, 10.0, 0.0];
        assert_eq!(
            cell_lines(&tr_bl, 4.0),
            vec![[(0.0, 0.4), (0.4, 0.0)], [(0.6, 1.0), (1.0, 0.6)]]
        );
        assert_eq!(
            cell_lines(&tr_bl, 6.0),
            vec![[(0.0, 0.6), (0.4, 1.0)], [(0.6, 0.0), (1.0, 0.4)]]
        );
    }

    /// The joiner this module used before edge keys: attach each segment to
    /// any line whose end point equals one of its points, then merge lines
    /// until no two share an end. Quadratic, but a straightforward reference
    /// for what the lines should be.
    fn chain_by_matching_points(segments: &[Segment]) -> Vec<(Vec<ContourPoint>, bool)> {
        let eq = |a: &ContourPoint, b: &ContourPoint| {
            (a.x - b.x).abs() < 1e-9 && (a.y - b.y).abs() < 1e-9
        };
        let mut chains: Vec<Vec<ContourPoint>> = Vec::new();
        for seg in segments {
            let attached = chains.iter_mut().any(|c| {
                let (first, last) = (c[0], c[c.len() - 1]);
                if eq(&last, &seg.p0) {
                    c.push(seg.p1);
                } else if eq(&last, &seg.p1) {
                    c.push(seg.p0);
                } else if eq(&first, &seg.p1) {
                    c.insert(0, seg.p0);
                } else if eq(&first, &seg.p0) {
                    c.insert(0, seg.p1);
                } else {
                    return false;
                }
                true
            });
            if !attached {
                chains.push(vec![seg.p0, seg.p1]);
            }
        }
        let mut merged = true;
        while merged {
            merged = false;
            'outer: for i in 0..chains.len() {
                for j in i + 1..chains.len() {
                    let (i_first, i_last) = (chains[i][0], chains[i][chains[i].len() - 1]);
                    let (j_first, j_last) = (chains[j][0], chains[j][chains[j].len() - 1]);
                    let mut taken = if eq(&i_last, &j_first)
                        || eq(&i_first, &j_last)
                        || eq(&i_last, &j_last)
                        || eq(&i_first, &j_first)
                    {
                        chains.remove(j)
                    } else {
                        continue;
                    };
                    if eq(&i_last, &j_first) {
                        taken.remove(0);
                        chains[i].append(&mut taken);
                    } else if eq(&i_first, &j_last) {
                        taken.pop();
                        taken.append(&mut chains[i]);
                        chains[i] = taken;
                    } else if eq(&i_last, &j_last) {
                        taken.reverse();
                        taken.remove(0);
                        chains[i].append(&mut taken);
                    } else {
                        taken.reverse();
                        taken.pop();
                        taken.append(&mut chains[i]);
                        chains[i] = taken;
                    }
                    merged = true;
                    break 'outer;
                }
            }
        }
        chains
            .into_iter()
            .map(|pts| {
                let closed = pts.len() >= 3 && eq(&pts[0], &pts[pts.len() - 1]);
                (pts, closed)
            })
            .collect()
    }

    /// A line as a canonical list of points, so two joiners' output can be
    /// compared whatever direction each walked a line in, and wherever each
    /// started a closed one.
    fn canonical(points: &[ContourPoint], closed: bool) -> Vec<(u64, u64)> {
        let mut pts: Vec<(u64, u64)> = points
            .iter()
            .map(|p| (p.x.to_bits(), p.y.to_bits()))
            .collect();
        if closed {
            pts.pop(); // the repeated first point
            let best = |v: &Vec<(u64, u64)>| -> Vec<(u64, u64)> {
                let k = (0..v.len()).min_by_key(|&i| v[i]).unwrap_or(0);
                v[k..].iter().chain(&v[..k]).copied().collect()
            };
            let forward = best(&pts);
            pts.reverse();
            let backward = best(&pts);
            forward.min(backward)
        } else {
            let mut reversed = pts.clone();
            reversed.reverse();
            pts.min(reversed)
        }
    }

    /// Every cell's segments at `level`, as `generate_contours` collects them.
    fn grid_segments(
        data: &[f64],
        width: usize,
        height: usize,
        level: f64,
        nodata: Option<f64>,
    ) -> Vec<Segment> {
        let mut segments = Vec::new();
        for row in 0..height - 1 {
            for col in 0..width - 1 {
                let corners = [
                    data[row * width + col],
                    data[row * width + col + 1],
                    data[(row + 1) * width + col],
                    data[(row + 1) * width + col + 1],
                ];
                if corners.iter().any(|&v| is_nodata(v, nodata)) {
                    continue;
                }
                let [tl, tr, bl, br] = corners;
                segments.extend(cell_segments(col, row, width, tl, tr, bl, br, level));
            }
        }
        segments
    }

    /// A deterministic grid of rolling terrain plus noise, with `holes`
    /// scattered nodata cells, so there are closed rings, lines that end at
    /// the border, and lines that end at a hole.
    fn bumpy_grid(width: usize, height: usize, seed: u64, holes: usize) -> Vec<f64> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut data: Vec<f64> = (0..width * height)
            .map(|i| {
                let (x, y) = ((i % width) as f64, (i / width) as f64);
                40.0 * (x / 5.0).sin() * (y / 7.0).cos()
                    + 15.0 * ((x + y) / 3.0).sin()
                    + 10.0 * next()
            })
            .collect();
        for _ in 0..holes {
            let i = (next() * (width * height) as f64) as usize % (width * height);
            data[i] = -9999.0;
        }
        data
    }

    #[test]
    fn test_chaining_by_edge_matches_chaining_by_point() {
        for (seed, holes) in [(1u64, 0usize), (2, 0), (3, 25), (4, 60)] {
            let (width, height) = (47, 39);
            let data = bumpy_grid(width, height, seed, holes);
            for level in [-30.0, -12.5, 0.0, 7.0, 21.0, 44.0] {
                let segments = grid_segments(&data, width, height, level, Some(-9999.0));
                assert!(
                    !segments.is_empty(),
                    "seed {seed} level {level} has no segments"
                );
                let mut by_edge: Vec<_> = chain_segments(&segments)
                    .iter()
                    .map(|(p, c)| (canonical(p, *c), *c))
                    .collect();
                let mut by_point: Vec<_> = chain_by_matching_points(&segments)
                    .iter()
                    .map(|(p, c)| (canonical(p, *c), *c))
                    .collect();
                by_edge.sort();
                by_point.sort();
                assert_eq!(
                    by_edge, by_point,
                    "seed {seed}, {holes} holes, level {level}"
                );
            }
        }
    }

    #[test]
    fn test_generate_contours_matches_per_level_reference() {
        // The whole pipeline against the way it used to run: every level
        // scanning every cell, then joining by matching points.
        let (width, height) = (53, 41);
        for (seed, holes) in [(5u64, 0usize), (6, 40)] {
            let data = bumpy_grid(width, height, seed, holes);
            for interval in [4.0, 9.5, 25.0] {
                let config = ContourConfig::new(interval)
                    .expect("valid config")
                    .with_base(1.25)
                    .with_nodata(-9999.0);
                let got = generate_contours(&data, width, height, &config).expect("contours");

                let valid: Vec<f64> = data
                    .iter()
                    .copied()
                    .filter(|v| !is_nodata(*v, config.nodata))
                    .collect();
                let (lo, hi) = valid
                    .iter()
                    .fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
                let mut want = Vec::new();
                for level in config.compute_levels(lo, hi) {
                    let segments = grid_segments(&data, width, height, level, config.nodata);
                    for (points, closed) in chain_by_matching_points(&segments) {
                        want.push((level.to_bits(), canonical(&points, closed), closed));
                    }
                }
                let mut got: Vec<_> = got
                    .iter()
                    .map(|c| {
                        (
                            c.level.to_bits(),
                            canonical(&c.points, c.is_closed),
                            c.is_closed,
                        )
                    })
                    .collect();
                got.sort();
                want.sort();
                assert!(!want.is_empty());
                assert_eq!(got, want, "seed {seed}, {holes} holes, interval {interval}");
            }
        }
    }

    #[test]
    fn test_chaining_uses_every_segment_once_when_values_equal_the_level() {
        // Integer elevations with a level they hit exactly: crossings land
        // on grid points, where several edges' points coincide. Lines may
        // touch there, but each segment belongs to exactly one line and
        // consecutive points of a line are the two ends of one segment.
        let (width, height) = (31, 23);
        let data: Vec<f64> = bumpy_grid(width, height, 9, 0)
            .iter()
            .map(|v| (v / 5.0).round() * 5.0)
            .collect();
        for level in [-20.0, 0.0, 10.0, 25.0] {
            let segments = grid_segments(&data, width, height, level, None);
            let lines = chain_segments(&segments);
            let used: usize = lines.iter().map(|(pts, _)| pts.len() - 1).sum();
            assert_eq!(
                used,
                segments.len(),
                "level {level}: every segment in exactly one line"
            );
            let mut pairs: Vec<_> = segments
                .iter()
                .map(|s| {
                    let (a, b) = (
                        (s.p0.x.to_bits(), s.p0.y.to_bits()),
                        (s.p1.x.to_bits(), s.p1.y.to_bits()),
                    );
                    (a.min(b), a.max(b))
                })
                .collect();
            let mut walked: Vec<_> = lines
                .iter()
                .flat_map(|(pts, _)| {
                    pts.windows(2).map(|w| {
                        let (a, b) = (
                            (w[0].x.to_bits(), w[0].y.to_bits()),
                            (w[1].x.to_bits(), w[1].y.to_bits()),
                        );
                        (a.min(b), a.max(b))
                    })
                })
                .collect();
            pairs.sort();
            walked.sort();
            assert_eq!(
                walked, pairs,
                "level {level}: lines are made of exactly the segments"
            );
            for (pts, closed) in &lines {
                if *closed {
                    assert_eq!(
                        pts.first(),
                        pts.last(),
                        "level {level}: a closed line ends where it starts"
                    );
                }
            }
        }
    }

    #[test]
    fn test_contours_md_worked_example() {
        // The worked example in CONTOURS.md, point for point.
        let data = [1.0, 2.0, 3.0, 2.0, 6.0, 4.0, 1.0, 3.0, 2.0];
        let config = ContourConfig::new(3.5).expect("valid config");
        let contours = generate_contours(&data, 3, 3, &config).expect("contours");
        assert_eq!(contours.len(), 1, "{contours:?}");
        let line = &contours[0];
        assert_eq!(line.level, 3.5);
        assert!(!line.is_closed);
        let want = [
            (2.0, 0.5),
            (1.0, 0.375),
            (0.375, 1.0),
            (1.0, 2.5 / 3.0 + 1.0),
            (2.0, 1.25),
        ];
        let got: Vec<(f64, f64)> = line.points.iter().map(|p| (p.x, p.y)).collect();
        let reversed: Vec<(f64, f64)> = got.iter().rev().copied().collect();
        let close = |a: &[(f64, f64)]| {
            a.iter()
                .zip(&want)
                .all(|(p, q)| (p.0 - q.0).abs() < 1e-12 && (p.1 - q.1).abs() < 1e-12)
        };
        assert!(
            got.len() == want.len() && (close(&got) || close(&reversed)),
            "{got:?}"
        );
    }

    #[test]
    fn test_contour_empty_grid() {
        let config = ContourConfig::new(10.0).expect("valid config");

        // Width 0
        let contours = generate_contours(&[], 0, 0, &config).expect("should succeed");
        assert!(contours.is_empty());

        // Height 0
        let contours = generate_contours(&[], 0, 5, &config).expect("should succeed");
        assert!(contours.is_empty());

        // Width 0, height > 0
        let contours = generate_contours(&[], 5, 0, &config).expect("should succeed");
        assert!(contours.is_empty());
    }

    #[test]
    fn test_contour_wrong_data_size() {
        let data = vec![1.0, 2.0, 3.0]; // 3 elements
        let config = ContourConfig::new(10.0).expect("valid config");

        // 2x2 expects 4 elements
        let result = generate_contours(&data, 2, 2, &config);
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_levels_with_base() {
        let config = ContourConfig::new(10.0)
            .expect("valid config")
            .with_base(5.0);
        let levels = config.compute_levels(0.0, 30.0);
        // base=5, interval=10, so levels: 5, 15, 25
        // But must be strictly inside (0, 30), so 5, 15, 25 all qualify
        assert!(levels.contains(&15.0));
        assert!(levels.contains(&25.0));
    }

    #[test]
    fn test_contour_config_nan_interval() {
        let result = ContourConfig::new(f64::NAN);
        assert!(result.is_err());
    }

    #[test]
    fn test_contour_config_inf_interval() {
        let result = ContourConfig::new(f64::INFINITY);
        assert!(result.is_err());
    }
}

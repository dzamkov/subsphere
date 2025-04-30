//! This test validates and maintains the `catalog.md` file.
use std::borrow::Cow;
use subsphere::prelude::*;

#[test]
fn test() {
    let mut cur = &*std::fs::read_to_string("catalog.md").unwrap();
    let mut act = String::new();

    // Process catalog line by line.
    let mut section = None;
    let mut is_valid = true;
    loop {
        // Keep track of section name
        if cur.starts_with("## [") {
            let end = cur.find("]").unwrap();
            section = Some(&cur[4..end]);
        }

        // Check for table header
        if cur.starts_with('|') {
            if let Some(section) = section {
                let (header_line, next) =
                    cur.split_at(cur.find('\n').expect("unexpected end of file") + 1);
                let header = parse_row(header_line).collect::<Box<[&str]>>();
                cur = next;
                let (_, next) = cur.split_at(cur.find('\n').expect("unexpected end of file") + 1);
                cur = next;
                let mut data = Vec::new();
                while cur.starts_with('|') {
                    let (row_line, next) = if let Some(line_end) = cur.find('\n') {
                        (&cur[..line_end + 1], &cur[line_end + 1..])
                    } else {
                        (cur, "")
                    };
                    let data_len = data.len();
                    data.extend(parse_row(row_line).map(Cow::Borrowed));
                    assert_eq!(
                        data.len() - data_len,
                        header.len(),
                        "row length mismatch in section `{}`",
                        section
                    );
                    cur = next;
                }
                check_table(section, &header, &mut data, &mut is_valid);

                // Rewrite the table
                let mut col_widths = header.iter().map(|s| s.len()).collect::<Box<_>>();
                for row in data.chunks_exact(header.len()) {
                    for (i, cell) in row.iter().enumerate() {
                        col_widths[i] = col_widths[i].max(cell.len());
                    }
                }
                write_row(&col_widths, header.iter().copied(), &mut act);
                act.push('|');
                for &col_width in &col_widths {
                    act.push(' ');
                    act.extend(std::iter::repeat_n('-', col_width));
                    act.push_str(" |");
                }
                act.push('\n');
                for row in data.chunks_exact(header.len()) {
                    write_row(&col_widths, row.iter().map(|s| s.as_ref()), &mut act);
                }
            }
        }

        // Advance to next line
        if let Some(line_end) = cur.find('\n') {
            act.push_str(&cur[..line_end + 1]);
            cur = &cur[line_end + 1..];
        } else {
            act.push_str(cur);
            break;
        }
    }

    // Overwrite the catalog if requested.
    if !is_valid {
        if std::env::var("OVERWRITE").is_ok() {
            std::fs::write("catalog.md", &act).unwrap();
        } else {
            panic!("'catalog.md' does not match output. Re-run with `OVERWRITE=1` to update it");
        }
    }
}

/// Splits a line of a table into its individual cells.
fn parse_row(line: &str) -> impl Iterator<Item = &str> {
    let mut items = line.split('|');
    items.next(); // Skip the prefix before the first '|'
    let mut items = items.map(str::trim);
    // Skip the last item
    let mut prev = items.next();
    std::iter::from_fn(move || {
        if let Some(cur) = items.next() {
            let res = prev.unwrap();
            prev = Some(cur);
            Some(res)
        } else {
            None
        }
    })
}

/// Writes a row of a table to the output string.
fn write_row<'a>(col_widths: &[usize], rows: impl Iterator<Item = &'a str>, out: &mut String) {
    out.push('|');
    for (i, cell) in rows.enumerate() {
        let width = col_widths[i];
        out.push_str(&format!(" {:width$} |", cell));
    }
    out.push('\n');
}

/// Validates the contents of a table.
///
/// If there is a discrepancy, this will set `is_valid` to `false` and update `rows`
/// to reflect the expected data.
fn check_table(section: &str, header: &[&str], data: &mut [Cow<str>], is_valid: &mut bool) {
    match section {
        "`TriSphere`" => {
            for row in data.chunks_exact_mut(header.len()) {
                let base = get_table_base("Base", header, row).expect("missing `Base` column");
                let b = std::num::NonZero::new(
                    get_table_str("B", header, row)
                        .expect("missing `B` column")
                        .parse::<u32>()
                        .unwrap(),
                )
                .unwrap();
                let c = get_table_str("C", header, row)
                    .expect("missing `C` column")
                    .parse::<u32>()
                    .unwrap();
                match get_table_str("Projector", header, row).expect("missing `Projector` column") {
                    "`Gnomonic`" => check_row(
                        subsphere::TriSphere::new(base, subsphere::proj::Gnomonic, b, c),
                        header,
                        row,
                        is_valid,
                    ),
                    "`Fuller`" => check_row(
                        subsphere::TriSphere::new(base, subsphere::proj::Fuller, b, c),
                        header,
                        row,
                        is_valid,
                    ),
                    str => panic!("unknown projector: {}", str),
                }
            }
        }
        "`HexSphere`" => {
            for row in data.chunks_exact_mut(header.len()) {
                let base = get_table_base("Base", header, row).expect("missing `Base` column");
                let b = std::num::NonZero::new(
                    get_table_str("B", header, row)
                        .expect("missing `B` column")
                        .parse::<u32>()
                        .unwrap(),
                )
                .unwrap();
                let c = get_table_str("C", header, row)
                    .expect("missing `C` column")
                    .parse::<u32>()
                    .unwrap();
                match get_table_str("Projector", header, row).expect("missing `Projector` column") {
                    "`Gnomonic`" => check_row(
                        subsphere::HexSphere::from_kis(subsphere::TriSphere::new(
                            base,
                            subsphere::proj::Gnomonic,
                            b,
                            c,
                        ))
                        .expect("invalid parameters"),
                        header,
                        row,
                        is_valid,
                    ),
                    "`Fuller`" => check_row(
                        subsphere::HexSphere::from_kis(subsphere::TriSphere::new(
                            base,
                            subsphere::proj::Fuller,
                            b,
                            c,
                        ))
                        .expect("invalid parameters"),
                        header,
                        row,
                        is_valid,
                    ),
                    str => panic!("unknown projector: {}", str),
                }
            }
        }
        _ => panic!("unknown section `{}` in catalog", section),
    }
}

/// Checks whether the values in a row correctly reflect the properties of the given sphere.
fn check_row<Sphere: subsphere::Sphere + Copy>(
    sphere: Sphere,
    header: &[&str],
    row: &mut [Cow<str>],
    is_valid: &mut bool,
) where
    Sphere::Face: std::fmt::Debug,
    Sphere::Vertex: std::fmt::Debug,
    Sphere::HalfEdge: std::fmt::Debug,
{
    // Validate basic consistency of the sphere
    subsphere::util::validate(sphere);

    // Compare measurements to expected values and update if necessary
    for (i, &col) in header.iter().enumerate() {
        let data = &mut row[i];
        match col {
            "# Faces" => {
                let exp = data.as_ref();
                let act = format!("{}", sphere.num_faces());
                if act != exp {
                    *is_valid = false;
                    *data = Cow::Owned(act);
                }
            }
            "# Verts" => {
                let exp = data.as_ref();
                let act = format!("{}", sphere.num_vertices());
                if act != exp {
                    *is_valid = false;
                    *data = Cow::Owned(act);
                }
            }
            "Area" => {
                let mut min_area = f64::MAX;
                let mut max_area = f64::MIN;

                // Only consider faces with maximum number of sides
                let num_sides = sphere.faces().map(|f| f.num_sides()).max().unwrap();
                for f in sphere.faces().filter(|f| f.num_sides() == num_sides) {
                    let area = f.area();
                    min_area = min_area.min(area);
                    max_area = max_area.max(area);
                }
                let exp = data.as_ref();
                let act = format!("{:.2}%", 100.0 * (max_area - min_area) / min_area);
                if act != exp {
                    *is_valid = false;
                    *data = Cow::Owned(act);
                }
            }
            "Length" => {
                let mut len_discrepancy: f64 = 0.0;
                for f in sphere.faces() {
                    let mut min_len = f64::MAX;
                    let mut max_len = f64::MIN;
                    for s in f.sides() {
                        let len = s.length();
                        min_len = min_len.min(len);
                        max_len = max_len.max(len);
                    }
                    len_discrepancy = len_discrepancy.max((max_len - min_len) / min_len);
                }
                let exp = data.as_ref();
                let act = format!("{:.2}%", 100.0 * len_discrepancy);
                if act != exp {
                    *is_valid = false;
                    *data = Cow::Owned(act);
                }
            }
            "Angle" => {
                let mut angle_discrepancy: f64 = 0.0;
                for f in sphere.faces() {
                    let mut min_angle = f64::MAX;
                    let mut max_angle = f64::MIN;
                    for s in f.sides() {
                        let angle = s.angle();
                        min_angle = min_angle.min(angle);
                        max_angle = max_angle.max(angle);
                    }
                    angle_discrepancy = angle_discrepancy.max(max_angle - min_angle);
                }
                let exp = data.as_ref();
                let act = format!("{:.2}°", angle_discrepancy * 180.0 / std::f64::consts::PI);
                if act != exp {
                    *is_valid = false;
                    *data = Cow::Owned(act);
                }
            }
            _ => (),
        }
    }
}

/// Gets the value of a table cell as a string.
fn get_table_str<'a>(column: &str, header: &[&str], row: &'a [Cow<str>]) -> Option<&'a str> {
    let index = header.iter().position(|&s| s == column)?;
    Some(row[index].as_ref())
}

/// Gets the value of a table cell as a [`subsphere::BaseTriSphere`].
fn get_table_base(
    column: &str,
    header: &[&str],
    row: &[Cow<str>],
) -> Option<subsphere::BaseTriSphere> {
    let value = get_table_str(column, header, row)?;
    match value {
        "`Icosa`" => Some(subsphere::BaseTriSphere::Icosa),
        "`Octo`" => Some(subsphere::BaseTriSphere::Octo),
        "`Tetra`" => Some(subsphere::BaseTriSphere::Tetra),
        _ => None,
    }
}

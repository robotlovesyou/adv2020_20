use std::collections::{HashMap, HashSet};
use std::{fmt, cmp};
use std::fmt::Formatter;
use std::rc::Rc;
use itertools::Itertools;

const TILE_SIZE: usize = 10;
const IMAGE_SIZE: usize = 12;
const BUFFER_SIZE: usize = 100;

fn main() {
    let tiles = read_tiles(include_str!("../input.txt").lines());
    let all_tiles = (&tiles[..]).all_variations();
    let part_1_solution = part_1(&all_tiles);
    println!("part 1: {}", part_1_solution);
}

struct Tile {
    id: u32,
    pixels: [u8; BUFFER_SIZE],
}

impl Tile {
    pub fn new(id: u32, pixels: [u8; BUFFER_SIZE]) -> Tile {
        Tile {
            id,
            pixels,
        }
    }

    pub fn connects_above(&self, other: &Self) -> bool {
        for x in 0..TILE_SIZE {
            if self.pixels[offset(x, TILE_SIZE-1)] != other.pixels[offset(x, 0)] {
                return false;
            }
        }
        true
    }

    pub fn connects_below(&self, other: &Self) -> bool {
        for x in 0..TILE_SIZE {
            if self.pixels[offset(x, 0)] != other.pixels[offset(x, TILE_SIZE-1)] {
                return false;
            }
        }
        true
    }

    pub fn connects_left_of(&self, other: &Self) -> bool {
        for y in 0..TILE_SIZE {
            if self.pixels[offset(TILE_SIZE-1, y)] != other.pixels[offset(0, y)] {
                return false;
            }
        }
        true
    }

    pub fn connects_right_of(&self, other: &Self) -> bool {
        for y in 0..TILE_SIZE {
            if self.pixels[offset(0, y)] != other.pixels[offset(TILE_SIZE-1, y)] {
                return false;
            }
        }
        true
    }


    pub fn flipped(&self) -> Rc<Self> {
        let mut new_pixels: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
        for y in 0..TILE_SIZE {
            for x in 0..TILE_SIZE {
                new_pixels[offset(TILE_SIZE-1-x, y)] = self.pixels[offset(x, y)];
            }
        }
        Rc::new(Tile::new(self.id, new_pixels))
    }


    pub fn rotated(&self) -> Rc<Self> {
        let mut rotated_pixels: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];

        for i in 0..TILE_SIZE {
            for j in 0..TILE_SIZE {
                rotated_pixels[offset(i, j)] = self.pixels[offset(TILE_SIZE-j-1,i)];
            }
        }
        Rc::new(Tile::new(self.id, rotated_pixels))
    }
}

trait Variable {
    type Item;
    fn variations(&self) -> Vec<Self::Item>;
}

impl Variable for Rc<Tile> {
    type Item = Rc<Tile>;
    fn variations(&self) -> Vec<Self::Item> {
        let mut variations = Vec::new();
        let mut current = self.clone();
        variations.push(current.clone());
        for _ in 0..3 {
            current = current.rotated();
            variations.push(current.clone());
        }

        current = self.flipped();
        variations.push(current.clone());
        for _ in 0..3 {
            current = current.rotated();
            variations.push(current.clone());
        }

        variations
    }
}

trait PictureOps {
    type Item;
    fn all_variations(&self) -> Vec<Self::Item>;
}

trait TilesOps {
    fn find_top_lefts(&self) -> Vec<Rc<Tile>>;
}

impl PictureOps for &[Rc<Tile>] {
    type Item = Rc<Tile>;
    fn all_variations(&self) -> Vec<Self::Item> {
        self.iter().map(|t| t.variations()).flatten().collect()
    }
}

impl TilesOps for &[Rc<Tile>] {
    fn find_top_lefts(&self) -> Vec<Rc<Tile>> {
        let mut top_lefts = Vec::new();

        'outer: for tile in self.iter() {
            for other in self.iter().filter(|t| t.id != tile.id) {
                if tile.connects_below(other) || tile.connects_right_of(other) {
                    continue 'outer;
                }
            }
            top_lefts.push(tile.clone());
        }
        top_lefts.iter().unique_by(|t| t.id).cloned().collect()
    }
}

impl fmt::Display for Tile {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut buffer = String::with_capacity(BUFFER_SIZE+TILE_SIZE);
        for y in 0..TILE_SIZE {
            for x in 0..TILE_SIZE {
                buffer += self.pixels[offset(x, y)].to_string().as_str();
            }
            buffer += "\n";
        }
        write!(f, "Tile{}\n{}", self.id, buffer)
    }
}

fn offset(x: usize, y: usize) -> usize {
    y * TILE_SIZE + x
}

fn read_tiles<'a>(lines: impl Iterator<Item=&'a str>) -> Vec<Rc<Tile>> {
    let mut peekable = lines.peekable();
    let mut tiles = Vec::new();
    let mut pixels: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    while peekable.peek().is_some() {
        let line = peekable.next().unwrap();
        let id = &line[5..line.len()-1].parse::<u32>().unwrap();

        for y in 0..TILE_SIZE {
            for (x, c) in peekable.next().unwrap().chars().enumerate() {
                if c == '#' {
                    pixels[offset(x, y)] = 1;
                }
            }

        }
        peekable.next().unwrap();
        tiles.push(Rc::new(Tile::new(*id, pixels)));
        pixels = [0; BUFFER_SIZE];
    }

    tiles
}

fn part_1(all_tiles: &[Rc<Tile>]) -> u64 {
    all_tiles.find_top_lefts().iter().map(|t| t.id as u64).product()
}


// #[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
// enum Edge {
//     Top,
//     Bottom,
//     Left,
//     Right
// }
//
// use Edge::*;
// use std::env::var;
//
// #[derive(Clone)]
// struct Tile {
//     id: u32,
//     pixels: Vec<Vec<char>>,
// }
//
// impl Tile {
//     pub fn new(id: u32, pixels: Vec<Vec<char>>) -> Tile {
//         Tile{id, pixels}
//     }
//
//     pub fn variations(self) -> Vec<Tile> {
//         let mut variations = Vec::new();
//         let mut current = self.flipped();
//         for _ in 0..3 {
//             let next = current.rotated();
//             variations.push(current);
//             current = next;
//         }
//         variations.push(current);
//
//         current = self;
//         for _ in 0..3 {
//             let next = current.rotated();
//             variations.push(current);
//             current = next;
//         }
//         variations.push(current);
//         variations
//     }
//
//     fn flipped(&self) -> Tile {
//         let mut flipped_pixels = self.pixels
//             .iter()
//             .map(|line| line.iter()
//                 .map(|c| *c)
//                 .rev()
//                 .collect())
//             .collect();
//
//         Tile::new(self.id, flipped_pixels)
//     }
//
//
//
//     fn rotated(&self) -> Tile {
//         let mut rotated_pixels = self.pixels.clone();
//         for i in 0..TILE_SIZE {
//             for j in 0..TILE_SIZE {
//                 rotated_pixels[i][j] = self.pixels[TILE_SIZE-j-1][i];
//             }
//         }
//         Tile::new(self.id, rotated_pixels)
//     }
//
//     fn top(&self) -> (String, Edge) {
//         (self.pixels[0].iter().collect(), Top)
//     }
//
//     fn bottom(&self) -> (String, Edge) {
//         (self.pixels[TILE_SIZE-1].iter().collect(), Bottom)
//     }
//
//     fn left(&self) -> (String, Edge) {
//         let mut left = String::new();
//         for i in 0..TILE_SIZE {
//             left.push(self.pixels[i][0]);
//         }
//         (left, Left)
//     }
//
//     fn right(&self) -> (String, Edge) {
//         let mut right = String::new();
//         for i in 0..TILE_SIZE {
//             right.push(self.pixels[i][TILE_SIZE-1]);
//         }
//         (right, Right)
//     }
// }
//
// impl fmt::Display for Tile {
//     fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
//         let mut out = format!("Tile: {}\n", self.id);
//         for i in 0..self.pixels.len() {
//             out += &format!("{}\n", self.pixels[i].iter().collect::<String>());
//         }
//         write!(f, "{}", out)
//     }
// }
//
// fn read_tiles<'a>(lines: impl Iterator<Item=&'a str>) -> Vec<Tile> {
//     let mut peekable = lines.peekable();
//     let mut tiles = Vec::new();
//     while peekable.peek().is_some() {
//         let line = peekable.next().unwrap();
//         let id = &line[5..line.len()-2].parse::<u32>().unwrap();
//         let mut pixels = Vec::new();
//         for _ in 0..TILE_SIZE {
//             let pixel_str = peekable.next().unwrap();
//             pixels.push(pixel_str.chars().collect());
//         }
//         peekable.next().unwrap();
//         tiles.push(Tile::new(*id, pixels))
//     }
//
//     tiles
// }
//
// fn main() {
//     let tiles = read_tiles(include_str!("../input.txt").lines());
//     let mut used: HashSet<u32> = HashSet::new();
//     let mut row_left = find_top_left(&tiles).expect("no top left");
//     // used.insert(top_left.id);
//     // println!("top left: {}", top_left);
//     //
//     //
//     // let mut left = top_left.clone();
//     // for _ in 0..11 {
//     //     let next_right = find_next_right(&tiles, &left, &used).expect("no next right");
//     //     used.insert(next_right.id);
//     //     println!("next_right: {}", next_right);
//     //     left = next_right;
//     // }
//     //
//     // let next_below = find_next_below(&tiles, &top_left, &used).expect("no next below");
//     // println!("{}", next_below);
//     let mut image: Vec<Vec<Tile>> = Vec::new();
//     for i in 0..IMAGE_SIZE {
//         let mut row = vec![row_left.clone()];
//         for _ in 0..11 {
//             let next_right = find_next_right(&tiles, &row[row.len() -1], &used).expect("no next right");
//             used.insert(next_right.id);
//             row.push(next_right);
//             println!("tile done");
//         }
//         if i < IMAGE_SIZE-1 {
//             row_left = find_next_below(&tiles, &row[0], &used).expect("no next below");
//         }
//         image.push(row);
//         println!("row done");
//     }
// }
//
// fn find_top_left(tiles: &[Tile]) -> Option<Tile> {
//     'outer: for tile in tiles.iter().map(|t| t.clone().variations()).flatten() {
//         let current = tile.clone();
//         let (top, _) = current.top();
//         let (left, _) = current.left();
//         let bottom_edge = (top, Bottom);
//         let right_edge = (left, Right);
//         for other in tiles.iter().filter(|t| t.id != current.id) {
//             for variant in other.clone().variations() {
//                 if bottom_edge == variant.bottom() {
//                     continue 'outer;
//                 }
//                 if right_edge == variant.right() {
//                     continue 'outer;
//                 }
//             }
//         }
//         return Some(current.clone());
//     }
//     None
// }
//
// fn find_next_right<'a>(tiles: &[Tile], left: &Tile, used: &HashSet<u32>) -> Option<Tile> {
//     let (right, _) = left.right();
//     let left_edge = (right, Left);
//     for other in tiles.iter().filter(|t| t.id != left.id && !used.contains(&t.id)) {
//         for variant in other.clone().variations() {
//             if left_edge == variant.left() {
//                 return Some(variant);
//             }
//         }
//     }
//     None
// }
//
// fn find_next_below<'a>(tiles: &[Tile], top: &Tile, used: &HashSet<u32>) -> Option<Tile> {
//     let (bottom, _) = top.bottom();
//     let top_edge = (bottom, Top);
//     for other in tiles.iter().filter(|t| t.id != top.id && !used.contains(&t.id)) {
//         for variant in other.clone().variations() {
//             if top_edge == variant.top() {
//                 return Some(variant);
//             }
//         }
//     }
//     None
// }

#[cfg(test)]
mod tests {
    use super::*;
    use indoc::indoc;

    const TEST_TILE: &str = indoc!("
        Tile 2311:
        ..##.#..#.
        ##..#.....
        #...##..#.
        ####.#...#
        ##.##.###.
        ##...#.###
        .#.#.#..##
        ..#....#..
        ###...#.#.
        ..###..###

        ");

    const TEST_TILE_EDGES: &str = indoc!("
        Tile 1234:
        .#.#.#.#.#
        ##..#.....
        ....##..##
        ####.#....
        .#.##.####
        ##...#.##.
        .#.#.#..##
        #.#....#..
        .##...#.##
        #.#.#.#.#.

        Tile 4567:
        #.#.#.#.#.
        .#..#.##.#
        ##.##.#...
        ..#.#.##.#
        #...#...#.
        ...##..###
        #..#.####.
        .#.####.##
        #.#..###..
        .#.#.#.#.#

        ");

    const TEST_QUERIES_AND_SOLUTIONS: &str = indoc!("
        Tile 2311:
        ..##.#..#.
        ##..#.....
        #...##..#.
        ####.#...#
        ##.##.###.
        ##...#.###
        .#.#.#..##
        ..#....#..
        ###...#.#.
        ..###..###

        Tile 1951:
        #.##...##.
        #.####...#
        .....#..##
        #...######
        .##.#....#
        .###.#####
        ###.##.##.
        .###....#.
        ..#.#..#.#
        #...##.#..

        Tile 1171:
        ####...##.
        #..##.#..#
        ##.#..#.#.
        .###.####.
        ..###.####
        .##....##.
        .#...####.
        #.##.####.
        ####..#...
        .....##...

        Tile 1427:
        ###.##.#..
        .#..#.##..
        .#.##.#..#
        #.#.#.##.#
        ....#...##
        ...##..##.
        ...#.#####
        .#.####.#.
        ..#..###.#
        ..##.#..#.

        Tile 1489:
        ##.#.#....
        ..##...#..
        .##..##...
        ..#...#...
        #####...#.
        #..#.#.#.#
        ...#.#.#..
        ##.#...##.
        ..##.##.##
        ###.##.#..

        Tile 2473:
        #....####.
        #..#.##...
        #.##..#...
        ######.#.#
        .#...#.#.#
        .#########
        .###.#..#.
        ########.#
        ##...##.#.
        ..###.#.#.

        Tile 2971:
        ..#.#....#
        #...###...
        #.#.###...
        ##.##..#..
        .#####..##
        .#..####.#
        #..#.#..#.
        ..####.###
        ..#.#.###.
        ...#.#.#.#

        Tile 2729:
        ...#.#.#.#
        ####.#....
        ..#.#.....
        ....#..#.#
        .##..##.#.
        .#.####...
        ####.#.#..
        ##.####...
        ##..#.##..
        #.##...##.

        Tile 3079:
        #.#.#####.
        .#..######
        ..#.......
        ######....
        ####.#..#.
        .#...#.##.
        #.#####.##
        ..#.###...
        ..#.......
        ..#.###...

    ");

    #[test]
    fn it_reads_a_single_tile() {
        assert_eq!(1, read_tiles(TEST_TILE.lines()).len());
    }

    #[test]
    fn it_checks_edge_connections() {
        let mut tiles = read_tiles(TEST_TILE_EDGES.lines());
        let tile_4567 = tiles.pop().unwrap();
        let tile_1234 = tiles.pop().unwrap();

        assert_eq!(4567, tile_4567.id);
        assert_eq!(1234, tile_1234.id);
        assert!(tile_1234.connects_above(&tile_4567));
        assert!(!tile_1234.connects_above(&tile_1234));

        assert!(tile_4567.connects_below(&tile_1234));
        assert!(!tile_4567.connects_below(&tile_4567));

        assert!(tile_1234.connects_left_of(&tile_4567));
        assert!(!tile_1234.connects_left_of(&tile_1234));

        assert!(tile_4567.connects_right_of(&tile_1234));
        assert!(!tile_4567.connects_right_of(&tile_4567));
    }

    #[test]
    fn a_flipped_tile_has_the_expected_pixels() {
        let tile = read_tiles(TEST_TILE.lines()).pop().unwrap();
        let flipped = tile.flipped();
        for y in 0..TILE_SIZE {
            for x in 0..TILE_SIZE {
                assert_eq!(tile.pixels[offset(x, y)], flipped.pixels[offset(TILE_SIZE-1-x, y)]);
            }
        }
    }

    #[test]
    fn a_rotation_of_a_tile_has_the_expected_properties() {
        let tile = read_tiles(TEST_TILE.lines()).pop().unwrap();
        let rotated = tile.rotated();
        assert_eq!(tile.id, rotated.id);

        let tile_left = (0..TILE_SIZE).map(|y| tile.pixels[offset(0, y)]).collect::<Vec<u8>>();
        assert_eq!(&tile_left[..], &rotated.pixels[BUFFER_SIZE-TILE_SIZE..]);

        let tile_right = (0..TILE_SIZE).map(|y| tile.pixels[offset(TILE_SIZE-1, y)]).collect::<Vec<u8>>();
        assert_eq!(&tile_right[..], &rotated.pixels[0..TILE_SIZE]);

        let rotated_left = (0..TILE_SIZE).map(|y| rotated.pixels[offset(0, y)]).collect::<Vec<u8>>();
        let tile_top = tile.pixels[0..TILE_SIZE].iter().rev().cloned().collect::<Vec<u8>>();
        assert_eq!(&rotated_left, &tile_top);

        let rotated_right = (0..TILE_SIZE).map(|y| rotated.pixels[offset(TILE_SIZE-1, y)]).collect::<Vec<u8>>();
        let tile_bottom = tile.pixels[BUFFER_SIZE-TILE_SIZE..].iter().rev().cloned().collect::<Vec<u8>>();
        assert_eq!(&rotated_right, &tile_bottom);
    }

    #[test]
    fn it_creates_the_expected_number_of_variations() {
        assert_eq!(8, read_tiles(TEST_TILE.lines()).pop().unwrap().variations().len());
    }

    #[test]
    fn it_finds_the_expected_top_lefts() {
        let all_tiles = (&read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines())[..]).all_variations();
        let top_lefts = (&all_tiles[..]).find_top_lefts();
        assert_eq!(4, top_lefts.len());
        assert!(top_lefts.iter().any(|t| t.id == 1951));
        assert!(top_lefts.iter().any(|t| t.id == 3079));
        assert!(top_lefts.iter().any(|t| t.id == 2971));
        assert!(top_lefts.iter().any(|t| t.id == 1171));
    }

    #[test]
    fn it_finds_the_expected_part_1_solution() {
        let all_tiles = (&read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines())[..]).all_variations();
        assert_eq!(20899048083289, part_1(&all_tiles));
    }
}
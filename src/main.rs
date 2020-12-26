use std::collections::{HashMap, HashSet};
use std::{fmt, cmp};
use std::fmt::Formatter;
use std::rc::Rc;

use im_rc::{HashSet as IMHashSet, Vector as IMVector, vector};
use itertools::{Itertools, all};

const TILE_SIZE: usize = 10;
const IMAGE_SIZE: usize = 12;
const IMAGE_PIXEL_SIZE: usize = IMAGE_SIZE * 8;
const BUFFER_SIZE: usize = 100;
const PICTURE_BUFFER_SIZE: usize = IMAGE_SIZE * IMAGE_SIZE * 8;

fn main() {
    let tiles = read_tiles(include_str!("../input.txt").lines());
    let all_tiles = (&tiles[..]).all_variations();
    let part_1_solution = part_1(&all_tiles);
    println!("part 1: {}", part_1_solution);

    let picture_vec = build_picture_vec(&all_tiles, IMAGE_SIZE);
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
    fn find_next_rights(&self, of: &Rc<Tile>, used: &IMHashSet<u32>) -> Vec<Rc<Tile>>;
    fn find_next_belows(&self, of: &Rc<Tile>, used: &IMHashSet<u32>) -> Vec<Rc<Tile>>;
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

    fn find_next_rights(&self, of: &Rc<Tile>, used: &IMHashSet<u32>) -> Vec<Rc<Tile>> {
        let mut next_rights = Vec::new();

        for tile in self.iter().filter(|t|!used.contains(&t.id) && t.id != of.id) {
            if tile.connects_right_of(of) {
                next_rights.push(tile.clone());
            }
        }

        next_rights
    }

    fn find_next_belows(&self, of: &Rc<Tile>, used: &IMHashSet<u32>) -> Vec<Rc<Tile>> {
        let mut next_belows = Vec::new();

        for tile in self.iter().filter(|t| !used.contains(&t.id) && t.id != of.id) {
            if tile.connects_below(of) {
                next_belows.push(tile.clone());
            }
        }

        next_belows
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

struct Image {
    pixels: Vec<u8>,
    size: usize,
}

fn image_offset(x: usize, y: usize, size: usize) -> usize {
    y * size + x
}

impl Image {
    fn new(tiles: &[Rc<Tile>], size: usize, tile_size: usize) -> Image {
        let mut pixels= vec![0u8; size * size];

        for (i, tile) in tiles.iter().enumerate() {
            let start_x = (i * 8) % (size);
            let start_y = (i / tile_size) * 8;
            println!("{},{}", start_x, start_y);
            for x in 1..TILE_SIZE-1 {
                for y in 1..TILE_SIZE-1 {
                    pixels[image_offset(start_x + x-1, start_y + y-1, size)] = tile.pixels[offset(x, y)];
                }
            }
        }

        Image {
            pixels,
            size
        }
    }
}

impl fmt::Display for Image {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut buffer = String::with_capacity(self.pixels.len());
        for y in 0..self.size {
            for x in 0..self.size {
                buffer += self.pixels[image_offset(x, y, self.size)].to_string().as_str();
            }
            buffer += "\n";
        }
        write!(f, "{}", buffer)
    }
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

fn build_picture_vec(all_tiles: &[Rc<Tile>], size: usize) -> Vec<Rc<Tile>> {
    let top_left = all_tiles.find_top_lefts().first().expect("no top lefts").clone();

    let result = picture_search(
        all_tiles,
        IMHashSet::new().update(top_left.id),
        vector![top_left.clone()],
        size
    ).expect("no result");

    let mut picture_vec = Vec::new();
    for tile in result.iter() {
        picture_vec.push(tile.clone());
    }
    picture_vec
}

fn picture_search(all: &[Rc<Tile>], used: IMHashSet<u32>, path: IMVector<Rc<Tile>>, size: usize) -> Option<IMVector<Rc<Tile>>> {
    if path.len() == size * size {
        return Some(path);
    }

    let candidates = match path.len() % size {
        0 => {
            let first_in_line = path.get(path.len() - size).expect("first in line");
            all.find_next_belows(first_in_line, &used.update(first_in_line.id))
        },
        _ => {
            let current = path.last().expect("no last");
            all.find_next_rights( current ,&used.update(current.id))
        },
    };

    for candidate in candidates {
        let mut new_path = path.clone();
        new_path.push_back(candidate.clone());
        if let Some(result) = picture_search(all, used.update(candidate.id), new_path, size) {
            return Some(result);
        }
    }
    None
}

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

    #[test]
    fn it_finds_the_expected_next_rights_of() {
        let tiles = read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines());

        // get tile 1951 in the correct orientation for tile 2311 in some rotation to be at its right
        let t1951_original = tiles.iter().find(|t|t.id == 1951).unwrap();
        let t1951 = t1951_original.flipped().rotated().rotated();

        let all_tiles = (&tiles[..]).all_variations();
        let right_of = (&all_tiles[..]).find_next_rights(&t1951, &IMHashSet::new());
        assert!(right_of.iter().find(|t| t.id == 2311).is_some());
    }

    #[test]
    fn it_finds_the_expected_next_below() {
        let tiles = read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines());

        // get tile 1951 in the correct orientation for tile 2311 in some rotation to be at its right
        let t1951_original = tiles.iter().find(|t|t.id == 1951).unwrap();
        let t1951 = t1951_original.flipped().rotated().rotated();

        let all_tiles = (&tiles[..]).all_variations();
        let below = (&all_tiles[..]).find_next_belows(&t1951, &IMHashSet::new());
        assert!(below.iter().find(|t| t.id == 2729).is_some());
    }

    #[test]
    fn it_builds_a_picture_from_the_test_input() {
        let all_tiles = (&read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines())[..]).all_variations();
        let picture_vec = build_picture_vec(&all_tiles, 3);
    }

    #[test]
    fn it_builds_an_image() {
        let all_tiles = (&read_tiles(TEST_QUERIES_AND_SOLUTIONS.lines())[..]).all_variations();
        let picture_vec = build_picture_vec(&all_tiles, 3);
        let image = Image::new(&picture_vec, 3 * 8, 3);
    }
}
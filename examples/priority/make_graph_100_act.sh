#!/bin/bash
# Creates /dev/shm/graph_100.act with 100 nodes and complex topology.
# Usage: ./make_graph_100_act.sh

cd "$(dirname "$0")/../.."

cat > /tmp/graph_100_gen.rs << 'EOF'
use mork::space::Space;
use std::fmt::Write;

fn main() {
    let mut s = Space::new();
    let mut data = String::new();

    // 1. The Ring (Nodes 0-49)
    for i in 0..50 {
        let next = (i + 1) % 50;
        writeln!(data, "(edge N{} N{})", i, next).unwrap();
    }

    // 2. The Hub (Node 50) - Connects to every 5th node in Ring
    for i in (0..50).step_by(5) {
        writeln!(data, "(edge N50 N{})", i).unwrap(); // Outgoing from Hub
        writeln!(data, "(edge N{} N50)", i + 2).unwrap(); // Incoming to Hub
    }

    // 3. The Tail (Nodes 51-70) - Linear chain attached to Ring at N25
    writeln!(data, "(edge N25 N51)").unwrap();
    for i in 51..70 {
        writeln!(data, "(edge N{} N{})", i, i + 1).unwrap();
    }

    // 4. The Island (Nodes 80-99) - Disconnected from Main Component
    // A small circle + random links
    for i in 80..99 {
        writeln!(data, "(edge N{} N{})", i, i + 1).unwrap();
    }
    writeln!(data, "(edge N99 N80)").unwrap(); // Close loop
    writeln!(data, "(edge N80 N85)").unwrap(); // Shortcut
    writeln!(data, "(edge N90 N95)").unwrap(); // Shortcut

    // 5. Incoming Connections (Island -> Main)
    // These make the Island visible IF we started there, but unreachable from N0.
    writeln!(data, "(edge N80 N49)").unwrap(); // Island Head -> Ring End
    writeln!(data, "(edge N99 N70)").unwrap(); // Island Tail -> Tail End

    println!("Generated {} edges.", data.lines().count());
    
    s.add_all_sexpr(data.as_bytes()).unwrap();
    s.backup_tree("/dev/shm/graph_100.act".to_string()).unwrap();
    println!("✓ Created /dev/shm/graph_100.act");
}
EOF

rustc --edition 2021 \
    -L target/release/deps \
    /tmp/graph_100_gen.rs \
    -o /tmp/graph_100_gen \
    --extern mork=target/release/libmork.rlib \
    --extern pathmap=$(ls target/release/deps/libpathmap-*.rlib | head -n 1)

/tmp/graph_100_gen

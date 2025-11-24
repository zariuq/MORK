#!/bin/bash
# Creates /dev/shm/graph_1.act for MORK
# Usage: ./make_graph_act.sh

# Move to MORK root relative to this script
cd "$(dirname "$0")/../.."

# Create temporary Rust source to generate the ACT
cat > /tmp/graph_act.rs << 'EOF'
use mork::space::Space;

fn main() {
    let mut s = Space::new();

    // Define the graph structure
    let data = r#"
(edge A B) (edge B C) (edge C A)
(edge A D) (edge D E) (edge E F) (edge F D)
(edge E F) (edge H G) (edge G A)
"#;

    s.add_all_sexpr(data.as_bytes()).unwrap();
    
    // Save to /dev/shm/ as required by MORK's default configuration
    s.backup_tree("/dev/shm/graph_1.act".to_string()).unwrap();
    println!("✓ Created /dev/shm/graph_1.act");
}
EOF

# Compile using the project's release libraries
rustc --edition 2021 \
    -L target/release/deps \
    /tmp/graph_act.rs \
    -o /tmp/graph_act \
    --extern mork=target/release/libmork.rlib \
    --extern pathmap=$(ls target/release/deps/libpathmap-*.rlib | head -n 1)

# Execute the generator
/tmp/graph_act

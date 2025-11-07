#!/bin/bash
# Creates /dev/shm/fruit_colors.act for ACT demo examples

cd "$(dirname "$0")/../.." # Go to MORK root

cat > /tmp/fruits_act.rs << 'EOF'
use mork::space::Space;

fn main() {
    let mut s = Space::new();

    let data = r#"
(color apple red)
(color banana yellow)
(color grape purple)
(color orange orange)
(color lemon yellow)
(color strawberry red)
"#;

    s.add_all_sexpr(data.as_bytes()).unwrap();
    s.backup_tree("/dev/shm/fruit_colors.act".to_string()).unwrap();
    println!("✓ Created /dev/shm/fruit_colors.act (6 color facts)");
}
EOF

rustc --edition 2021 -L target/release/deps /tmp/fruits_act.rs \
  -o /tmp/fruits_act \
  --extern mork=target/release/libmork.rlib \
  --extern pathmap=target/release/deps/libpathmap-*.rlib

/tmp/fruits_act

#!/bin/bash
# Creates /dev/shm/math_facts.act for ACT demo examples

cd "$(dirname "$0")/../.." # Go to MORK root

cat > /tmp/math_act.rs << 'EOF'
use mork::space::Space;

fn main() {
    let mut s = Space::new();

    let data = r#"
(add 0 0 0)
(add 0 1 1)
(add 0 2 2)
(add 1 0 1)
(add 1 1 2)
(add 1 2 3)
(add 1 3 4)
(add 2 0 2)
(add 2 1 3)
(add 2 2 4)
(add 2 3 5)
(add 3 0 3)
(add 3 1 4)
(add 3 2 5)
(add 3 3 6)
"#;

    s.add_all_sexpr(data.as_bytes()).unwrap();
    s.backup_tree("/dev/shm/math_facts.act".to_string()).unwrap();
    println!("✓ Created /dev/shm/math_facts.act (15 addition facts)");
}
EOF

rustc --edition 2021 -L target/release/deps /tmp/math_act.rs \
  -o /tmp/math_act \
  --extern mork=target/release/libmork.rlib \
  --extern pathmap=target/release/deps/libpathmap-*.rlib

/tmp/math_act

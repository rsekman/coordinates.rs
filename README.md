# Coordinates.rs

## Getting started

Add the following to your `cargo.toml` under `[dependencies]`

```toml
coordinates = "0.4.0"
```

if you want additional features — like serializing and deserializing — your `[dependencies]` line will look more like this
```toml
coordinates = { version = "0.4.1", features = ["serde"] }
```

In a file import the coordinate system you want, or all of them through 
`coordinates::prelude::*`

```rust
use coordinates::two_dimensional::Vector2;
```
If you want extra traits, such as magnitude or dot products you will also need 
to include the following

```rust
use coordinates::traits::*;
```

And finally initialize a variable

```rust
let var = Vector2 {
  x: 0.0,
  y: 1.0,
}
```

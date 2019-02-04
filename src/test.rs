#[cfg(test)]
mod cip_internal;
mod tests {
    let current = GridParameters::new(
        f: vec![0.1],
        u: vec![0.1],
        g: vec![0.1],
        dx: 1.0,
        dt: 1.0
    )
    #[test]
    fn it_works(current: current) {
        assert_eq!(current.f, vec![0.1]);
    }
}

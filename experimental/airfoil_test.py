import pterasoftware as ps

# Individual tests (sanity visualization):
# bad:
# okay: e337, fx62k131, mh150, strand
test_airfoil = ps.geometry.airfoil.Airfoil(
    name="strand",
)

test_airfoil.draw()

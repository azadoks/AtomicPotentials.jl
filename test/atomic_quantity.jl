# TODO: test `ht`
# TODO:    for analytical quantities -> just change the type parameter
# TODO:    for numerical quantities -> perform the Hankel transform and store
# TODO:        the new mesh & function

# TODO: test `iht`
# TODO:    ditto ^

# TODO: test `interpolate_onto`
# TODO:     for analytical quantities -> evaluate the quantity on the given mesh, return
# TODO:         the corresponding type that supports numerical quantities
# TODO:     for numerical quantities -> evaluate the quantitie's interpolator on the given
# TODO:         mesh, return a new struct of the same type with the new mesh &
# TODO:         function-on-the-mesh
# TODO:         NOTE: this could (maybe more intuitively) be done using a `Base.convert`
# TODO:             dispatch so that interpolate_onto doesn't perform two conceptually
# TODO:             different things, i.e. (interpolating) evaluating on a (new) mesh

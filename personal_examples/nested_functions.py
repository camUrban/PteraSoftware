new_values = []


def outer_function(outer_args):
    print("In outer function.")
    print("Outer args from outer function:" + str(outer_args))

    def inner_function(inner_args):
        print("In inner function.")
        print("Outer args from inner function:" + str(outer_args))
        print("Inner args from inner function: " + str(inner_args))

        new_values.append("Whoa")

    inner_function(inner_args=outer_args + ["Hello"])
    print(new_values)


outer_function(["My", "Name", "Is", "Cameron"])

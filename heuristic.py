def next_step(coordinate, direction, movement):
    if direction == "U":
        if movement == "f":
            next_coordinate = (coordinate[0], coordinate[1] + 1)
            new_direction = "U"
        elif movement == "r":
            next_coordinate = (coordinate[0]+1, coordinate[1])
            new_direction = "R"
        elif movement == "l":
            next_coordinate = (coordinate[0]-1, coordinate[1])
            new_direction = "L"
    if direction == "R":
        if movement == "f":
            next_coordinate = (coordinate[0]+1, coordinate[1])
            new_direction = "R"
        elif movement == "r":
            next_coordinate = (coordinate[0], coordinate[1]-1)
            new_direction = "D"
        elif movement == "l":
            next_coordinate = (coordinate[0], coordinate[1]+1)
            new_direction = "U"
    if direction == "D":
        if movement == "f":
            next_coordinate = (coordinate[0], coordinate[1]-1)
            new_direction = "D"
        elif movement == "r":
            next_coordinate = (coordinate[0]-1, coordinate[1])
            new_direction = "L"
        elif movement == "l":
            next_coordinate = (coordinate[0]+1, coordinate[1])
            new_direction = "R"
    if direction == "L":
        if movement == "f":
            next_coordinate = (coordinate[0]-1, coordinate[1])
            new_direction = "L"
        elif movement == "r":
            next_coordinate = (coordinate[0], coordinate[1]+1)
            new_direction = "U"
        elif movement == "l":
            next_coordinate = (coordinate[0], coordinate[1]-1)
            new_direction = "D"
    return next_coordinate, new_direction


def find_closest_path(current_path, index, possible_options):
    """
    Function to find the closest path of an index that contains a 1 to a possible pairing position
    Takes into account the current path to avoid overlapping and to determine the number of steps that need to be taken

    Returns the path from the last point in the sequence to the possible pairing point
    If no such path exist, return None
    """


def find_hp(input_string):
    """
    Function to get a hp-folding in 2D using a greedy heuristic of pairing the next 1 to the closest possible position
    This function will assume a path allowing 3 movements (f,l,r)
    """
    #Get all the even and odd positions with ones
    even_1 = [i  for i, value in enumerate(input_string) if i%2 == 0 and value == "h"]
    odd_1 = [i  for i, value in enumerate(input_string) if i%2 == 1 and value == "h"]
    all_1 = [i  for i, value in enumerate(input_string) if value == "h"]

    #Possible pairing positions for even and odd 1's
    pairing_odd = []
    pairing_even = []

    #Initialization of path, starts by convention in coordinate (0,0) and points to the right
    last_point = (0,0)
    current_direction = "R"
    current_index = 0

    #path with coordinates of the folding
    path = [last_point]

    #Sequence of movements to get the fold
    movements = ""

    while current_index < len(input_string):
        if len(pairing_odd) == 0 and len(pairing_even) == 0:
            movements += "f"
            last_point, current_direction = next_step(last_point, current_direction, "f")
            path.append(last_point)
            if current_index in all_1:
                current_index = all_1.pop(0)
                possible_neighbors = [(last_point[0]-1, last_point[1]),
                                        (last_point[0]+1, last_point[1]),
                                        (last_point[0], last_point[1]+1),
                                        (last_point[0], last_point[1]-1)]
                if current_index%2 == 0:
                    even_1.pop(0)
                    neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_odd]
                    pairing_odd += neighbors
                else:
                    odd_1.pop(0)
                    neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_even]
                    pairing_even += neighbors
            current_index += 1
            print(pairing_even, pairing_odd)

    return movements

#input_string = "hhppppphhppphppphp"
input_string = "pppppppppppppppphh"
print(find_hp(input_string))
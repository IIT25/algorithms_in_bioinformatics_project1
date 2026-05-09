from collections import deque
from time import time

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


def find_closest_path(current_path, current_coord, current_dir, possible_options, steps_allowed, priority):
    """
    Function to find the closest path of an index that contains a 1 to a possible pairing position
    Takes into account the current path to avoid overlapping and to determine the number of steps that need to be taken

    Returns the path from the last point in the sequence to the possible pairing point
    If no such path exist, return None
    """
    if not possible_options or steps_allowed <= 0:
        return None

    # Queue stores: (coord, direction, movement_sequence, visited_set)
    # Using a set for current_path lookup is O(1)
    initial_visited = set(current_path)
    queue = deque([(current_coord, current_dir, "", initial_visited)])

    while queue:
        coord, direction, moves, visited = queue.pop()
        # If we've reached the required length
        if len(moves) == steps_allowed:
            if coord in possible_options:
                if (coord[0]-1, coord[1]) not in visited or (coord[0]+1, coord[1]) not in visited or \
                    (coord[0], coord[1]+1) not in visited or (coord[0], coord[1]-1) not in visited:                                
                    return moves
            continue # Can't go further than steps_allowed

        # Try all 3 relative movements
        for m in priority:
            nxt_coord, nxt_dir = next_step(coord, direction, m)

            # Check for self-intersection
            if nxt_coord not in visited:
                new_visited = visited | {nxt_coord} # Set union
                queue.append((nxt_coord, nxt_dir, moves + m, new_visited))

    return None


def find_hp(input_string, priority):
    """
    Function to get a hp-folding in 2D using a greedy heuristic of pairing the next 1 to the closest possible position
    This function will assume a path allowing 3 movements (f,l,r)
    """
    #Get all the even and odd positions with ones
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

    if 0 in all_1:
        current_index = all_1.pop(0)
        possible_neighbors = [(last_point[0]-1, last_point[1]),
                                (last_point[0]+1, last_point[1]),
                                (last_point[0], last_point[1]+1),
                                (last_point[0], last_point[1]-1)]
        if current_index%2 == 0:
            neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_odd]
            
            pairing_odd += neighbors
        else:
            neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_even]
            pairing_even += neighbors


    while current_index < len(input_string)-1:
        if len(pairing_odd) == 0 and len(pairing_even) == 0:
            movements += "f"
            last_point, current_direction = next_step(last_point, current_direction, "f")
            current_index += 1
            path.append(last_point)
            if current_index in all_1:
                current_index = all_1.pop(0)
                possible_neighbors = [(last_point[0]-1, last_point[1]),
                                        (last_point[0]+1, last_point[1]),
                                        (last_point[0], last_point[1]+1),
                                        (last_point[0], last_point[1]-1)]
                if current_index%2 == 0:
                    neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_odd]
                    pairing_odd += neighbors
                else:
                    neighbors = [possible_neighbor for possible_neighbor in possible_neighbors if possible_neighbor not in path and possible_neighbor not in pairing_even]
                    pairing_even += neighbors
            
        else:
            """Condition to look forward, need to find the closest point in the string after the current index
                that can be located in a point of pairing even or pairing odd respecting the conditions.
                To do this, can look into odd_1 and comapre with pairing odd, or look into even_1 and compare with pairing_even"""
            found_path = None

            # Look forward at upcoming 'h' indices
            for h_idx in all_1:
                # Parity Rule: Even index H pairs with Odd index neighbors and vice versa
                targets = pairing_odd if h_idx % 2 != 0 else pairing_even
                steps_needed = h_idx - current_index
                if steps_needed == 1:
                    last_neighbors = [(last_point[0]-1, last_point[1]),
                                        (last_point[0]+1, last_point[1]),
                                        (last_point[0], last_point[1]+1),
                                        (last_point[0], last_point[1]-1)]
                    while targets and targets[-1] in last_neighbors:
                        targets = targets[:-1]
                # Search for a path to any valid pairing neighbor
                path_str = find_closest_path(path, last_point, current_direction, targets, steps_needed, priority)
                if path_str:
                    found_path = path_str
                    break

            if found_path:
                # Execute the movements found
                for move in found_path:
                    movements += move
                    last_point, current_direction = next_step(last_point, current_direction, move)
                    path.append(last_point)
                    current_index += 1
                    
                    if current_index in all_1:
                        # Logic to update pairing_even/odd for the newly placed H
                        possible_neighbors = [(last_point[0]-1, last_point[1]),
                                            (last_point[0]+1, last_point[1]),
                                            (last_point[0], last_point[1]+1),
                                            (last_point[0], last_point[1]-1)]

                        # If target_h_index is even, its neighbors are candidates for odd H's
                        if current_index % 2 == 0:
                            new_neighbors = [n for n in possible_neighbors if n not in path and n not in pairing_odd]
                            pairing_odd += new_neighbors
                        else:
                            new_neighbors = [n for n in possible_neighbors if n not in path and n not in pairing_even]
                            pairing_even += new_neighbors

                        # Clean up all_1 list (remove the H we just used and any we skipped)
                        all_1 = [i for i in all_1 if i > current_index]
                        break

            else:
                # Fallback: No path to a pairing, just step to the next available space
                for movement in ["f", "l", "r"]:

                    point, direction = next_step(last_point, current_direction, movement)
                    if point not in path:
                        last_point = point
                        current_direction = direction
                        movements += movement
                        path.append(last_point)
                        current_index += 1
                        break
                    return movements
    return movements

def find_possible_hp(input_string):
    priorities = [["f", "l", "r"], ["r", "f", "l"], ["r", "l", "f"]]
    results = []
    for priority in priorities:
        result = find_hp(input_string, priority)
        if len(result) == len(input_string)-1:
            results.append(result)
    return results

inputs = ["hhppppphhppphppphp", "hphphhhppphhhhpphh", "phpphphhhphhphhhhh", "hphpphhphpphphhpphph", "hhhpphphphpphphphpph", "hhpphpphpphpphpphpphpphh",
          "pphpphhpppphhpppphhpppphh", "ppphhpphhppppphhhhhhhpphhpppphhpphpp", "pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh",
          "hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh", "pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp",
          "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh", "hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph",
          "pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh",
          "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"]
for input_string in inputs:
    print(len(input_string))
    t0 = time()
    results = find_possible_hp(input_string)
    t1 = time()
    print(t1-t0)
    with open("heuristic results.txt", "a") as f:
        f.write(str(input_string) + " got the following results: "+ str(results) +" in time " + str(t1-t0) + "\n")
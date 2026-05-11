from time import time

def get_even_odd_lists(string: str) -> tuple:
    """
    Args:
        string (str): sequence of letters h and p

    Returns:
        tuple: two lists with indexes for even's and odd's positions of letter h
    """
    even = []
    odd = []
    for i in range(len(string)):
        char = string[i]
        if char == "h" or char == "H":
            if i % 2 == 0:
                even.append(i)
            else:
                odd.append(i)
    return even, odd

def find_best_folding_point(string: str, even: list, odd: list) -> tuple:
    """
    Args:
        string (str): sequence of letters h and p
        even (list): even's position of letter h
        odd (list): odd's position of letter h

    Returns:
        tuple: 
            best size of maximum matching
            best folding point
            best match type (even left/odd right (EL/OR) or odd left/even right (OL/ER))
    """
    best_size = -1
    best_point = -1
    best_match_type = ""
    repeated_count = 1
    for index in range(len(string)):
        even_left = sum(1 for i in even if i < index)
        even_right = sum(1 for i in even if i > index+1)
        odd_left = sum(1 for i in odd if i < index)
        odd_right = sum(1 for i in odd if i > index+1)
    
        match1 = min(even_left, odd_right)
        match2 = min(odd_left, even_right)
        
        current_max = max(match1, match2)
        if current_max == 0:
            continue
        if current_max == match1:
            left_first = max(i for i in even if i < index)
            right_first = min(i for i in odd if i > index+1)
        else:
            left_first = max(i for i in odd if i < index)
            right_first = min(i for i in even if i > index+1)
        
        best = (right_first+left_first)//2
        
        if current_max > best_size:
            repeated_count = 1
            best_size = current_max
            best_point = best
            if match1 == current_max:
                best_match_type = "EL/OR"
            else:
                best_match_type = "OL/ER"
        elif current_max == best_size:
            repeated_count+= 1
    if repeated_count%2 == 0:
        best_point += repeated_count//2-1
    
    return best_size, best_point, best_match_type

def generate_fold(string: str, even: list, odd: list, best_size: int, best_point: int, best_match_type: str):
    moves = ""
    
    if best_match_type == "EL/OR":
        left_block = [even_idx for even_idx in even if even_idx < best_point]
        right_block = [odd_idx for odd_idx in odd if odd_idx > best_point+1]
    else:        
        left_block = [odd_idx for odd_idx in odd if odd_idx < best_point]
        right_block = [even_idx for even_idx in even if even_idx > best_point+1] 

    #Left block
    for current_index in range(left_block[0]):
        moves += "f"
    
    for block_index in range(len(left_block) - 1):
        current_point = left_block[block_index]
        next_pairing_point = left_block[block_index+1]
        dist = next_pairing_point - current_point
        if dist == 2:
            moves += "ff"
        else:
            moves += "l" + "f"*(dist//2-2) + "rr" + "f"*(dist//2-2) + "l"
    
    steps_to_fold = 0
    for current_index in range(left_block[-1], best_point):
        moves += "f"
        steps_to_fold += 1

    # 180-degree turn
    moves += "rr"
    
    #Right block
    for j in range(best_point+2, right_block[0]):
        moves += "f" 
    
    for block_index in range(len(right_block) - 1):
        current_point = right_block[block_index]
        next_pairing_point = right_block[block_index+1]
        dist = next_pairing_point - current_point
        if dist == 2:
            moves += "ff"
        else:
            moves += "l" + "f"*(dist//2-2) + "rr" + "f"*(dist//2-2) + "l"
    
    for current_index in range(right_block[-1], len(string)-1):
        moves += "f"

    return moves

inputs = ["hhppppphhppphppphp", "hphphhhppphhhhpphh", "phpphphhhphhphhhhh", "hphpphhphpphphhpphph", "hhhpphphphpphphphpph", "hhpphpphpphpphpphpphpphh",
          "pphpphhpppphhpppphhpppphh", "ppphhpphhppppphhhhhhhpphhpppphhpphpp", "pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh",
          "hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh", "pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp",
          "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh", "hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph",
          "pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh",
          "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"]

for string in inputs:
    t0 = time()
    even, odd = get_even_odd_lists(string)
    best_size, best_point, best_match_type = find_best_folding_point(string, even, odd)
    result = generate_fold(string, even, odd, best_size, best_point, best_match_type)
    t1 = time()
    with open("approximation results.txt", "a") as f:
        f.write(str(string) + " got the following result: "+ str(result) +" in time " + str(t1-t0) + "\n")
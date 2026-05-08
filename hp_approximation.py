

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
    for index in range(len(string)):
        even_left = sum(1 for i in even if i <= index)
        even_right = len(even) - even_left
        odd_left = sum(1 for i in odd if i <= index)
        odd_right = len(odd) - odd_left
    
        match1 = min(even_left, odd_right)
        match2 = min(odd_left, even_right)
        
        current_max = max(match1, match2)
        
        if current_max > best_size:
            best_size = current_max
            best_point = index
            if match1 == current_max:
                best_match_type = "EL/OR"
            else:
                best_match_type = "OL/ER"
    return best_size, best_point, best_match_type
    
string = "hhhhhhhhhhhhphphpp"
even, odd = (get_even_odd_lists(string))
print(find_best_folding_point(string, even, odd))

### If your loop finds two different indices that give the same best_size, the algorithm usually prefers the one that sits on a Polar (P) bead or a long sequence of zeros, as this provides more "room" to make the turn on the lattice.
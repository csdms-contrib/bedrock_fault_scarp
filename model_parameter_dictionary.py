#! /usr/env/python

class ModelParameterDictionary:
    """Reads model parameters from an input file to a dictionary
       and provides functions for the user to look up particular parameters
       by key name."""
    #------------------------------------------------------------
    # The __init__ method just creates an empty dictionary
    # caled param_dict
    #------------------------------------------------------------
    def __init__(self):
    	self.param_dict = {}
    #------------------------------------------------------------
    # read_from_file: given the name of an input file, read
    # parameters. The format of the input file should be like:
    #
    #   # A comment line
    #   SOME_KEY: this the key for some parameter
    #   1.234
    #
    #   In other words, the rules are:
    #    - Comments are preceded by hash characters
    #    - Each parameter has two consecutive lines, one for the
    #      key and one for the value
    #    - The key must be followed by a space, colon, or eol
    #    - The parameter can be numeric or text
    #
    #------------------------------------------------------------
    def read_from_file(self, file_name):
        try:
            input_file = open( file_name )
        except IOError:
            print 'Unable to open', file_name
            raise
        stripped_line_list = []
        for line in input_file:
            line = line.strip()   # strip leading spaces
            if len(line)>0 and line[0] != '#':
                stripped_line_list.append( line )
        iskey = True
        for line in stripped_line_list:
            if iskey:
                # Strip out everything after the first space or colon
                first_colon = line.find(':')
                if first_colon == -1: first_colon = len(line)
                first_space = line.find(' ')
                if first_space == -1: first_space = len(line)
                last_key_char = min( first_colon, first_space )
                last_key = line[0:last_key_char]  # remember the key
                iskey = False
            else:
                self.param_dict[last_key] = line
                iskey = True
        input_file.close()
    #------------------------------------------------------------
    def handle_missing_key( self, key ):
        print 'I could not find the key "', key, '" in the input file.'
        raise KeyError
    #------------------------------------------------------------
    def handle_value_error( self, key, expected_type ):
        print 'In reading the value for "' + key + '",',
        print 'I expected', expected_type + '.'
        print 'Instead I found "' + str( self.param_dict[key] ) + '".'
        raise ValueError
    #------------------------------------------------------------
    # read_int:
    # Use this like: i = read_int( 'MY_INT' )
    # An error is generated if MY_INT isn't in the dictionary or
    # if its value is not an integer.
    #------------------------------------------------------------
    def read_int(self, key):
        try: my_value = self.param_dict[key]
        except KeyError: self.handle_missing_key( key )
        try: my_int = int( my_value )
        except ValueError: self.handle_value_error( key, 'an integer' )
        return my_int
    #------------------------------------------------------------
    # read_float:
    # Use this like: x = read_float( 'MY_FLOAT' )
    # An error is generated if MY_FLOAT isn't in the dictionary or
    # if its value is not a number.
    #------------------------------------------------------------
    def read_float(self, key):
        try: my_value = self.param_dict[key]
        except KeyError: self.handle_missing_key( key )
        try: my_float = float( my_value )
        except ValueError:
            self.handle_value_error( key, 'a floating-point number' )
        return my_float
    #------------------------------------------------------------
    # read_string:
    # Use this like: s = read_string( 'MY_STRING' )
    # An error is generated if MY_STRING isn't in the dictionary.
    #------------------------------------------------------------
    def read_string(self, key):
        try: my_value = self.param_dict[key]
        except KeyError: self.handle_missing_key( key )
        return str( my_value )
    #------------------------------------------------------------
    # read_bool:
    # Use this like: b = read_bool( 'MY_BOOL' )
    # An error is generated if MY_BOOL isn't 0, 1, True or False
    #------------------------------------------------------------
    def read_bool(self, key):
        try: my_value = self.param_dict[key]
        except KeyError: self.handle_missing_key( key )
        if my_value=='True' or my_value=='1' or my_value==1:
            return True
        elif my_value=='False' or my_value=='0' or my_value==0:
            return False
        else:
            self.handle_value_error( key, 
                              'a boolean (0, 1, True or False)' )
    #------------------------------------------------------------
    # read_int_cmdline: reads an integer from the command line
    # Use this like: i = read_int_cmdline( 'MY_INT' )
    # An error is generated if MY_INT is not an integer.
    #------------------------------------------------------------
    def read_int_cmdline(self, key):
        my_value = input( key + ': ' )
        self.param_dict[key] = my_value
        if isinstance( my_value, int ) != True:
            self.handle_value_error( key, 'an integer' )
        return my_value
    #------------------------------------------------------------
    # read_float_cmdline: reads a float from the command line
    # Use this like: f = read_int_cmdline( 'MY_FLOAT' )
    # An error is generated if MY_FLOAT is not a float.
    #------------------------------------------------------------
    def read_float_cmdline(self, key):
        my_value = input( key + ': ' )
        self.param_dict[key] = my_value
        try: my_float = float( my_value )
        except ValueError:
            self.handle_value_error( key, 'a floating-point number' )
        return my_float
    #------------------------------------------------------------
    # read_string_cmdline: reads a float from the command line
    # Use this like: s = read_string_cmdline( 'MY_STRING' )
    #------------------------------------------------------------
    def read_string_cmdline(self, key):
        my_str = raw_input( key + ': ' )
        self.param_dict[key] = my_str
        return my_str
    #------------------------------------------------------------
    # read_bool_cmdline:
    # Use this like: b = read_bool_cmdline( 'MY_BOOL' )
    # An error is generated if MY_BOOL isn't 0, 1, True or False
    #------------------------------------------------------------
    def read_bool_cmdline(self, key):
        my_value = raw_input( key + ': ' )
        if my_value=='True' or my_value=='1' or my_value==1:
            return True
        elif my_value=='False' or my_value=='0' or my_value==0:
            return False
        else:
            self.handle_value_error( key, 
                              'a boolean (0, 1, True or False)' )
        
                
            

            
            
        






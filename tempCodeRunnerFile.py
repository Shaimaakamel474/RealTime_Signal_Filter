    def list_of_lists_to_complex(self,complex_data):
        if isinstance(complex_data[0], list):
            complex_numbers = [complex(r, imag) for r, imag in complex_data]
        else:
            complex_numbers = complex(complex_data[0], complex_data[1])
        return complex_numbers
    
    def from_str_to_float(self,selected_mode_str):
        minus_sign_count = selected_mode_str.count('-')
        plus_sign_count = selected_mode_str.count('+')
        parts = re.split(pattern=r"([-()+j])", string=selected_mode_str)
        numeric_parts = [part for part in parts if part and part not in {'-', '(', '+', ')', 'j'}]
        if minus_sign_count == 2:
            numeric_floats = [-float(part) for part in numeric_parts]
        else:
            if minus_sign_count == 1 and plus_sign_count:
                numeric_floats = [float(part) for part in numeric_parts]
                numeric_floats[0]=  -abs(numeric_floats[0])
            else:
                numeric_floats = [float(part) for part in numeric_parts]
                numeric_floats[-1]=  -abs(numeric_floats[-1]) 
        
        if minus_sign_count == 0:
            numeric_floats = [float(part) for part in numeric_parts]

        if len(numeric_floats) == 1:
            numeric_floats.insert(0, 0.0)
        return numeric_floats
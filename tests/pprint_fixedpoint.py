import gdb

# Pretty printer for FixedPoint class
class FixedPointPrinter:
    "Pretty Printer for FixedPoint<T, FractionalBits>"

    def __init__(self, val):
        self.val = val
        self.template_params = self.extract_template_params()

    # Strip typedefs and aliases to get the underlying FixedPoint<T, FractionalBits> type
    def extract_template_params(self):
        # Strip away any typedefs or aliases to get to the underlying type
        real_type = self.val.type.strip_typedefs()

        self.type_str = str(real_type)
        # Example self.type_str: 'FixedPoint<int64_t, 60>'
        start = self.type_str.find('<') + 1
        end = self.type_str.find('>')
        if start != -1 and end != -1:
            params = self.type_str[start:end].split(', ')
            if len(params) == 3:
                T_type = gdb.lookup_type(params[0])
                LONG_T_type = gdb.lookup_type(params[0])
                fractional_bits = int(params[2])
                return T_type, LONG_T_type, fractional_bits
        return None, None

    # This function calculates the floating-point value
    def to_float(self):
        T_type, LONG_T_type, fractional_bits = self.template_params
        if T_type and LONG_T_type and fractional_bits is not None:
            # Extract the raw fixed-point value (T value)
            fixed_value = int(self.val['value'])
            # Perform fixed-point to float conversion: value / (1 << FractionalBits)
            return float(fixed_value) / (1 << fractional_bits)
        return None

    # This function returns the floating-point representation as a string
    def to_string(self):
        float_value = self.to_float()
        if float_value is not None:
            return "FixedPoint<float> = {:.6f}".format(float_value)
        return f"FixedPoint<unresolved>: {self.type_str}"

# Function to register the pretty printer
def build_pretty_printers():
    pp = gdb.printing.RegexpCollectionPrettyPrinter("my_pretty_printers")
    # Match FixedPoint with two template parameters
    pp.add_printer('FixedPoint', '^FixedPoint<.*>$', FixedPointPrinter)
    gdb.printing.register_pretty_printer(gdb.current_objfile(), pp)

# Call the registration function automatically when the script is loaded
build_pretty_printers()
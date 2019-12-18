abstract type AbstractXCCalculator
end

# Use internal XC
struct XCCalculator <: AbstractXCCalculator
end

# Use Libxc
struct LibxcXCCalculator <: AbstractXCCalculator
end

# These types should contain xcfunc field

uT = u"yr"
uL = u"m"
uM = u"g"
uQ = uM*uL^(-2)  # unit of pools, i.e. quantities: g/m^2

ustrip_num(num::Num) = Unitful.ustrip(ModelingToolkit.get_unit(num), ModelingToolkit.value(num))

#ustrip(::)
# Unitful.ustrip(ModelingToolkit.get_unit(k),v)
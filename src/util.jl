"""
    smoothstep(x, x_step, dx, a=zero(x), b=one(x))

smooth step function: 
- = a for x <= x_step - dx
- = b for x >= x_step + dx
- = (a-b)/2 for x == x_step
and smooth in between.
"""
function smoothstep(x, x_step, dx, a=zero(x), b=one(x))
    @assert dx > 0
    edge0 = x_step - dx
    x <= edge0 && return(a)
    edge1 = x_step + dx
    x >= edge1 && return(b)
    y0 = _smoothstep(x,edge0,edge1)
    y0 * (b-a) + a
end

function _smoothstep(x, edge0, edge1) 
    #https://en.wikipedia.org/wiki/Smoothstep
    # Scale, bias and saturate x to 0..1 range
    # assume that 
    @assert edge0 < x < edge1
    #x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
    xs = (x - edge0) / (edge1 - edge0)
    # Evaluate polynomial
    xs * xs * (3 - 2 * xs)
  end
  

  


  
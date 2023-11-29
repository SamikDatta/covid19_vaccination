function CIString = makeCI(x, scale, form)

xScaled = x/scale;

CIString = sprintf( [form, ' [', form, ', ', form, ']'], xScaled([2, 1, 3]) );



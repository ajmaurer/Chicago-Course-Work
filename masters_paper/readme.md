This contains work related to my masters paper, where I worked on a small extension the knockoff filter method developed by Emmanuel Candès and my advisor Rina Foygel Barber.

Please note, the code is still in a rough form and its likely a few of the options and specifications are broken. Please feel free to contact me if you're interested and have questions or issues.

Abstract:
Variable selection for regression is a key problem in applied statistics. The knockoff filter method provides one method of variable selection for linear regression. It relies on generating 'knockoff' features, which replicate the correlation structure of the original variable; when the full path of LASSO regression is fit, the point at which a null variable and its knockoff have nonzero coefficients will be exchangeable. However, for other GLMs, the method breaks down. I will provide an alternative method of randomly generating knockoffs for binary variables which will will satisfy the original correlation condition in expectation and over improved performance for other GLMs.

1. In the "gene_expression_count" excel file
2. I first found the mean values of the Control and treatment values and in some values we get zeroes to counter this the following strategies are used

Fold Chain is calculated by 
= avgerage of Treatment/ average of control values

if the value is undeterminable (0/0) we change it to "Nan"

later to get logFC we apply log base 2 to our Fold Change values

here also we face a similar issue of indertiminable terms(log (Nan), log (0))
so we have used the below strategies to encounter the issue:
1. if the values are "Nan" or "0" we used "Nan" as the result for the problem


later while plotting the error bar graphs we consider the Nan points value equal to zero.

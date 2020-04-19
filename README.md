# ThermalConduction2D
A model for solving 2-D heat conduction through multiple materials with user defined IVP's.

These are the initial notes I made for what I would like to include:

-Material represented by 2D array

--	each cell represents dX.dY (done)

---		calculated at start (array)

--	properties unique to each cell

---		properties can be calculated for each cell

--	can be imported / designed by user

---		use excel spreadsheet + macros, or image tool such as paint.net to "paint" the materials and IVP

---   write a GUI instead (TKinter?)


-Initial Value Problem and Calculations

-- 	calculate for time T (done)

--- 	use variety of solving methods

--  allow for properties to be re-calculated on each iteration instead of assuming constant


-Return results as a visual representation

-- 	Heatmap can show each property specified by user/pre-determined (done)

--  Sub-plot of start and end condition 

--  Make use of animation to show change of heat spread/properties with time?


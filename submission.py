How to match fields from two lists and further filter based upon the values in subsequent fields?



I am attempting to get the pos and alt strings from file1 to match up with what is in 
file2, fairly simple. However, file2 has values in the 17th split element/column to the 
last element/column (340th) which contains string such as 1/1:1.2.2:51:12 which 
I also want to filter for. 

I want to extract the rows from file2 that contain/match the pos and alt from file1. 
Thereafter, I want to further filter the matched results that only contain certain 
values in the 17th split element/column onwards. But to do so the values would have to 
be split by ":" so I can filter for split[0] = "1/1" and split[2] > 50. The problem is 
I have no idea how to do this.

I imagine I will have to iterate over these and split but I am not sure how to do this 
as the code is presently in a loop and the values I want to filter are in columns not rows.

Any advice would be greatly appreciated, I have sat with this problem since Friday and 
have yet to find a solution.

    import os,itertools,re
    file1 = open("file1.txt","r")
    file2 = open("file2.txt","r")

    matched = []

    for (x),(y) in itertools.product(file2,file1):
        if not x.startswith("#"):
                cells_y = i.split("\t")
                pos_y = cells[0]
                alt_y = cells[3]

                cells_x = x.split("\t")
                pos_x = cells_x[0]+":"+cells_x[1]
                alt_x = cells_x[4]

                if pos_y in pos_x and alt_y in alt_x:
                        matched.append(x)

    for z in matched:
        cells_z = z.split("\t")
        if cells_z[16:len(cells_z)]:


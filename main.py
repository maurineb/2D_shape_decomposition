import os
import subprocess


# OPTIONS
# Path of the image
# image_path = "rat.png"
# image_path = "test_star.png"
# image_path = "test_loop_outside_main.png"
# image_path = "test_different_loops.png"
# image_path = "other_test_loop_outside_main.png"
image_path = "lizzard-9.png"
# image_path = "apple-18.png"
# Value of epsilon for the skeleton
epsilon = "1.0"
# Show edf, wedf, et, st like matlab
matlab = "false"





# Create a folder for the results
commande = "mkdir results"
subprocess.run(commande, shell=True)
os.chdir(".\\results")
commande = "mkdir skeleton"
subprocess.run(commande, shell=True)
commande = "mkdir bma"
subprocess.run(commande, shell=True)
commande = "mkdir decomposition"
subprocess.run(commande, shell=True)
commande = "mkdir values"
subprocess.run(commande, shell=True)

# Go to the folder for the generation of the skeleton
commande = "..\\compact-skel-2d"
os.chdir(commande)

# Compute the skeleton
commande = "cargo run --release -- --imgin ..\\" + image_path + " --imgout ..\\results\\skeleton\\skeleton.png" + " --epsilon " + epsilon + " --nodout ..\\results\\skeleton\\nod.txt --edgout ..\\results\\skeleton\\edg.txt --delout ..\\results\\skeleton\\del.txt --bndout ..\\results\\skeleton\\bnd.txt" 
subprocess.run(commande, shell=True)

# Go to the folder for the decomposition
commande = "..\\decomposition\\bma_decomposition"
os.chdir(commande)

# Compute the decomposition
commande = "cargo run --release -- --folder_skel ..\\..\\results\\skeleton --folder_results ..\\..\\results"
subprocess.run(commande, shell=True)


from PIL.Image import *

nom_image = "../ressources/test_star.png"
i = open(nom_image)

(largeur, hauteur)= i.size
for x in range(largeur):
    for y in range(hauteur):
        (rouge,vert,bleu, alpha) = i.getpixel((x,y))
        if rouge+vert+bleu < 127.5:
            i.putpixel((x,y),(0,0,0,alpha))
        else :
            i.putpixel((x,y),(255,255,255,alpha))

i.save("../ressources/test_star.png")


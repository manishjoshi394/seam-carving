# Seam-Carving
**What's this ?** Get the literature reference here. https://en.wikipedia.org/wiki/Seam_carving

**See the API documentation here: http://manishjoshi394.github.io/seam-carving**

**Seam carving** (or **liquid rescaling**) is an algorithm for content-aware image resizing, developed by **Shai Avidan**, of **Mitsubishi Electric Research Laboratories (MERL)**, and **Ariel Shamir**, of the **Interdisciplinary Center** and **MERL**. It functions by establishing a number of seams (paths of least importance) in an image and automatically removes seams to reduce image size or inserts seams to extend it.  
It resizes the Image keeping in mind the objects in the image.
Here's an example,

- ***Original Image to be made narrower :***

![Original Chameleon](/chameleon.png)


- ***Scaling is Undesirable because chameleon is distorted :***

![Scaled Chameleon](/docs/img/rescale_chameleon.png)


- ***Cropping is Undesirable because either chameleon or branch and the credit text is removed :***

![Cropped Chameleon](/docs/img/crop_chameleon.jpg)

- ***This is what Seam-Carving achieves :***

![Carved Chameleon](/chameleon_resized.png)

Seam Carving takes into account the most noticeable parts of the image to us, preserves such parts; removes the less important pixels from the image. 

This repository contains Java based API for Seam-Carving. The API allows the client to edit images using the methods provided by it. 
## How to use the API ?
- Download a copy of the repository.
- Copy the *SeamCarver.java*, dependencies folder and LICENSE file to your Project's working directory. 
- Now you can import *SeamCarver.java* and include it directly in your program.
- You should be able to create a `SeamCarver` object by providing Image location as a constructor argument.
- You can use the methods `findVerticalSeam()` and `removeVerticalSeam(int[] seam)` etc. provided by `SeamCarver` class to edit the Image as required.
- Use the `picture()` method to get a `Picture` object corresponding to resized image. You can store it using a `Picture` class reference.
- Use the method `save(File file)` in thus stored `Picture` object to store the Picture to a disk.
- You may use the method `show()` in the `Picture` object to show the Image onscreen.
- You can use sample demo clients `ResizeDemo.java` and `SCUtility.java` to test how the images get resized. Read the source code headers of these files for compilation instructions. 

**Note: `Picture` class used here is a part of the Educational [algs4 library](https://github.com/kevin-wayne/algs4). Read It's API Documentation [here](https://algs4.cs.princeton.edu/code/javadoc/edu/princeton/cs/algs4/Picture.html)**

**Work is Going on..**
**ReadMe under construction**
## API reference: http://manishjoshi394.github.io/seam-carving

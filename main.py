from PyQt5 import QtCore, QtWidgets, uic
import sys 
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import matplotlib.pyplot as plt
import pydicom
import pydicom.data
from PIL import Image
import os
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import magic
import numpy as np
import math 
from collections import Counter
import random
from random import randint
import sys
from PyQt5.QtCore import Qt


def quicksort(array):
    if len(array) < 2:
        return array
    low, same, high = [], [], []
    pivot = array[randint(0, len(array) - 1)]
    for item in array:
        if item < pivot:
            low.append(item)
        elif item == pivot:
            same.append(item)
        elif item > pivot:
            high.append(item)
    return quicksort(low) + same + quicksort(high)

def rescale_intensities(image_arr, width, height):
        # rescale the pixel intensities to be between 0 and 255
        for i in range(height):
            for j in range(width):
                if image_arr[i][j] < 0:
                    image_arr[i][j] = 0
                elif image_arr[i][j] > 255:
                    image_arr[i][j] = 255
        return image_arr


def zero_pad(kernel_size, img_height, img_width, image_array): # add zero padding to the image
    kernel_pad = kernel_size // 2 # thickness of padding "frame"
    padded_image = np.zeros((img_height + (kernel_size - 1), img_width + (kernel_size - 1))) # setting array of the size of the padded image
    # starting at the position where the first image pixel is, which is going to be the kernel center
    # ending at the image dimension + the amount padded
    for i in range(kernel_pad, img_height + kernel_pad):
        for j in range(kernel_pad, img_width + kernel_pad):
            padded_image[i][j] = image_array[i - kernel_pad][j - kernel_pad] # inserting image data within the frame of the padding
    return padded_image
    

class Ui(QtWidgets.QMainWindow): # class to load .ui in constructor
    def __init__(self):
        super(Ui, self).__init__() # call the inherited classes __init__ method
        uic.loadUi('image_viewer.ui', self) # load the .ui file
        self.setWindowTitle("Image Viewer")
        self.show() # Show the GUI
        # for image viewer (DICOm and normal images)
        self.scene = QGraphicsScene()
        self.image_view.setScene(self.scene) #set a graphics scene in QGraphicsView
        # for original image in rotation tab
        self.scene1 = QGraphicsScene()
        self.original_view.setScene(self.scene1) #set a graphics scene in QGraphicsView
        # for rotated image
        self.scene2 = QGraphicsScene()
        self.edited_view.setScene(self.scene2) #set a graphics scene in QGraphicsView
        # for histogram of opened image
        self.scene3 = QGraphicsScene()
        self.orig_hist_view.setScene(self.scene3) #set a graphics scene in QGraphicsView
        #for histogram of equalized image
        self.scene4 = QGraphicsScene()
        self.equal_hist_view.setScene(self.scene4) #set a graphics scene in QGraphicsView
        # for magnitude of FFT
        self.scenem = QGraphicsScene()
        self.magnitude_view.setScene(self.scenem) #set a graphics scene in QGraphicsView
        # for logged magnitude of FFT
        self.sceneml = QGraphicsScene()
        self.magnitude_log_view.setScene(self.sceneml) #set a graphics scene in QGraphicsView
        # for phase of FFT
        self.sceneph = QGraphicsScene()
        self.phase_view.setScene(self.sceneph) #set a graphics scene in QGraphicsView
        # for logged phase of FFT
        self.scenephl = QGraphicsScene()
        self.phase_log_view.setScene(self.scenephl) #set a graphics scene in QGraphicsView

        self.histogram_flag = 0

        self.action_open.triggered.connect(self.browse_files) # call browse_files when Open is clicked in menubar
        self.zoom_button.clicked.connect(self.zooming)
        self.plot_T()
        self.nearest_button.clicked.connect(self.rotate_nearest)
        self.bilinear_button.clicked.connect(self.rotate_bilinear)
        self.shear_button.clicked.connect(self.shear)
        self.action_open.triggered.connect(self.show_image_tab4)
        self.action_open.triggered.connect(self.show_image_tab5)
        self.unsharp_button.clicked.connect(self.highboost_enhancement)
        self.shear_button.clicked.connect(self.shear) 
        self.orig_histogram.clicked.connect(self.original_histogram)
        self.equal_histogram.clicked.connect(self.equalized_histogram)
        self.unsharp_button.clicked.connect(self.highboost_enhancement)
        self.season_button.clicked.connect(self.seasoning)
        self.unseason_button.clicked.connect(self.median_filter)
        self.clear_median_button.clicked.connect(self.clear_median_label)
        self.clear_highboost_button.clicked.connect(self.clear_higboost_label)
        self.action_open.triggered.connect(self.show_image_tab7)
        self.action_open.triggered.connect(self.fourier_transform)
                   
    def check_tab(self):
        return self.tab_widget.currentIndex() 
        
    def clear_tab_1(self):
        self.width_line.clear() 
        self.height_line.clear() 
        self.size_line.clear()
        self.color_line.clear()
        self.depth_line.clear()
        self.modality_line.clear()
        self.name_line.clear()
        self.age_line.clear()
        self.part_line.clear()
        self.tab4_check = 0

    def browse_files(self): #browse device to open any file
        try:
            self.file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File') #qt function to open file browser, no specified type so it opens all files
            # the , _  next to file_name allows the first parameter to receive the path of the image file as a str 
            # rather than a tuple which was not working with the Qpixmap function and displayed the error: TypeError: QPixmap(): argument 1 has unexpected type 'tuple'
            # path = file_name[0] #variable to store opened file path
            splitter = os.path.split(self.file_path) # splitting file path from name to be used in dicom function
            self.file_destination = splitter[0]
            self.file_name = splitter[1]
            self.status_bar.showMessage(self.file_path)  # shows file path in status bar
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                self.dicom() 
            else: # call the normal image and data functions (for jpg, bmp, etc.)
                self.show_image()
                self.obtain_data()
        except:
            # clearing past entries
            self.clear_tab_1()
            self.scene.clear()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('Please choose a file to open')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def show_image(self): # shows any common image type (jpg, bmp, png, etc.)
        if self.check_tab() == 0:
            try:
                self.scene.clear() # clears the scene of the past image to show new image
                self.image = QImage(self.file_path) #places the selected file as the image QImage class allows direct access to the pixel data to represent an image
                pic = QGraphicsPixmapItem() # QGraphicsPixmapItem class to create a pixmap that can added to a QGraphicsScene, constructing an item
                pic.setPixmap(QPixmap.fromImage(self.image)) # QImage conversion into a QPixmap using fromImage() and setting it as the Pixmap
                self.image_view.fitInView(pic, QtCore.Qt.KeepAspectRatio) # make image fit in the Graphics View while keeping the aspevct ratio
                self.scene.addItem(pic) # add the image item to the created scene 
            except:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("ERROR")
                msg.setInformativeText('An error has occured')
                msg.setWindowTitle("Error")
                msg.exec_()
    
    def obtain_data(self): # gets and displays the data related to the opened image
        #getting the data

        try:
            image = Image.open(self.file_path) # read image data to find image.shape which has an output (height, width, number of channels) and if it is one channel (height, width)
            # determine the number of channels with the understanding of the output from image.shape
            np_img = np.asarray(image)
            if len(np_img.shape) == 2:
                self.num_channels = 1
            else:
                self.num_channels = np_img.shape[-1]

            img = Image.open(self.file_path) # open image to get min and max pixel values
            min = np.amin(image)
            max = np.amax(image)
            depth = math.ceil(math.log((int(max) - int(min) + 1), 2)) 
            depth_channels = self.num_channels * depth # calculate depth using min, max, and number of channels
            width, height = img.size # get image width and height
            if depth_channels == 0:
                depth_channels = 1
                depth = 1
            size_bits = width * height * depth # calculate size in bits

            # determine image color (RGB, Greyscale, Binary) from min, max, and number of channels
            if self.num_channels == 3:
                self.color = "Color (RGB)"
            elif self.num_channels == 4:
                self.color = "Color (RGBA)"
            elif self.num_channels == 1 and depth == 8:
                self.color = "Greyscale"
            elif self.num_channels == 2 and depth == 8:
                self.color = "Greyscale"
            elif self.num_channels == 1 and depth == 1:
                self.color = "Binary"
            else:
                self.color = "Other"
          
            # clearing past entries
            self.clear_tab_1()

            #displaying the data in GUI
            self.width_line.setText(str(width)) 
            self.height_line.setText(str(height)) 
            self.size_line.setText(str(size_bits)) 
            self.color_line.setText(self.color)
            self.depth_line.setText(str(depth_channels))

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def dicom(self): # displays the dicom image and reads and displays the image data
        
        try:
            # reading the .dcm file data
            filename = pydicom.data.data_manager.get_files(self.file_destination, self.file_name)[0] # get dicom file using path and name seperately
            ds = pydicom.dcmread(filename) # reading the data in the .dcm file
            self.scene.clear()
            # show the dicom image in the GUI       
            figure = Figure() # using matplotlib's Figure which holds all plot elements
            ax = figure.gca() # get the current axes
            # remove appearance of axes to show the image
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.imshow(ds.pixel_array, cmap=plt.cm.bone) # take pixel data from dicom file and display the data as an image in the bone color map
            canvas = FigureCanvas(figure) # pass the figure to Figure canvas which handles all the details of talking to user interface toolkits (pyqt in this case)
            self.scene.addWidget(canvas) # add a widget to the scene in QGraphicsView and use canvas to interface with matplotlib
            self.image_view.setScene(self.scene) # set the scene in the ui's graphics view
        
            # clearing past displayed data
            self.clear_tab_1()

            # displaying required data with condition in case an attribute is not in a .dcm file
            if hasattr(ds, 'Columns'):
                self.width_line.setText(str(ds.Columns))
            else: 
                self.width_line.setText("DOES NOT EXIST")

            if hasattr(ds, 'Rows'):
                self.height_line.setText(str(ds.Rows))
            else: 
                self.height_line.setText("DOES NOT EXIST")
            
            if hasattr(ds, 'BitsAllocated'):
                self.depth_line.setText(str(ds.BitsAllocated))
            else: 
                self.depth_line.setText("DOES NOT EXIST")
            
            file_size = ds.Columns * ds.Rows * ds.BitsAllocated
            self.size_line.setText(str(file_size))
            
            if hasattr(ds, 'BitsAllocated'):
                self.depth_line.setText(str(ds.BitsAllocated))
            else: 
                self.depth_line.setText("DOES NOT EXIST")
            
            if hasattr(ds, 'PhotometricInterpretation'):
                self.color_line.setText(str(ds.PhotometricInterpretation))
            else: 
                self.color_line.setText("DOES NOT EXIST")
            
            if hasattr(ds, 'Modality'):
                self.modality_line.setText(str(ds.Modality))
            else: 
                self.modality_line.setText("DOES NOT EXIST")

            if hasattr(ds, 'PatientName'):
                self.name_line.setText(str(ds.PatientName))
            else: 
                self.name_line.setText("DOES NOT EXIST") 

            if hasattr(ds, 'PatientAge'):
                self.age_line.setText(str(ds.PatientAge))
            else: 
                self.age_line.setText("DOES NOT EXIST") 

            if hasattr(ds, 'StudyDescription'):
                self.part_line.setText(str(ds.StudyDescription))
            else: 
                self.part_line.setText("DOES NOT EXIST")  
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def nearest_neighbor(self, image_array, height, width, new_height, new_width, factor):
        new_arr = np.zeros((new_height, new_width)) # 2D array of zeros with size same as zoomed image
        # replace the zeros with the correct pixel data 
        for i in range(new_height):
            for j in range(new_width):
                if i/factor > (height - 1) or j/factor > (width -1):
                    new_arr[i,j] =  image_array[int(i/factor), int(j/factor)]
                else:
                    new_arr[i,j] =  image_array[round(i/factor), round(j/factor)]

        # display image
        im = Image.fromarray(np.uint8(new_arr))
        qimg = im.toqpixmap() # create qpixmap to display image in Qlabel
        self.neighbor_label.clear() # clear past images
        self.neighbor_label.setPixmap(qimg) # set the created pixmap in the label
    
    def bilinear_interpolation(self, image_arr, y, x):
        height = image_arr.shape[0]
        width = image_arr.shape[1]

        x1 = max(min(math.floor(x), width - 1), 0) # max(sth,0) to avoid negative indeces, min(x,width - 1) & min(x,height - 1) to stay within original image bounds
        y1 = max(min(math.floor(y), height - 1), 0) # floor to get x1 and y1
        x2 = max(min(math.ceil(x), width - 1), 0) # ceil to get x2 and y2
        y2 = max(min(math.ceil(y), height - 1), 0)
        
        # getting the indices of 4 pixels in the original image that surround the pixel I am interpolating in the new image 
        # to calculate the new pixel data
        a = float(image_arr[y1, x1]) 
        b = float(image_arr[y2, x1])
        c = float(image_arr[y1, x2]) 
        d = float(image_arr[y2, x2]) 

        # distance between pixel I am interpolating and a
        dx = x - x1 
        dy = y - y1

        #calculating new pixel value after interpolation by multiplying old values at a,b,c,d by the distance using linear interpolation principles
        new_pixel = a * (1 - dx) * (1 - dy) # 1 - dx because here dx will be zero, same for dy
        new_pixel += b * dy * (1 - dx) # 1 - dx because here dx will be zero
        new_pixel += c * dx * (1 - dy) # 1 - dy because here dy will be zero
        new_pixel += d * dx * dy 
        return round(new_pixel)

    def resize_image(self, image_arr, new_height, new_width, factor):
        new_image = np.zeros((new_height, new_width)) # 2D array of zeros with size same as zoomed image

        orig_height = image_arr.shape[0]
        orig_width = image_arr.shape[1]

        # Find center column and center row
        x_orig_center = (orig_width-1) / 2
        y_orig_center = (orig_height-1) / 2

        # Find center of resized image
        x_scaled_center = (new_width-1) / 2
        y_scaled_center = (new_height-1) / 2

        # Compute the scale in both axes
        scale = 1 / factor 

        # fixed point resizing like in computer graphics we're scaling from new_size to old_size to pass it to the bilinear interpolation function
        # Computing x_scaled and y_scaled by:
        # Subtracting x_scaled_center and y_scaled_center from x and y which are the loop indices
        # After subtraction, (0, 0) is the "new center"
        # Scale the "zero centered" coordinates by multiplying the scale
        # Convert the "scaled zero centered" coordinates to "top left (0, 0)" by adding x_orig_center and y_orig_center
        for y in range(new_height):
            for x in range(new_width):
                # loop on the 2D array and perform bilinear interpolation on each pixel
                x_scaled = (x - x_scaled_center) * scale + x_orig_center
                y_scaled = (y - y_scaled_center) * scale + y_orig_center
                new_image[y, x] = self.bilinear_interpolation(image_arr, y_scaled, x_scaled)              
        return new_image

    def zooming(self):
        try:     
            factor = float(self.zooming_line.text())
            if factor <= 0: # check for negative or zero factor an disolay error message
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("ERROR")
                msg.setInformativeText('Factor must be greater than zero')
                msg.setWindowTitle("Error")
                msg.exec_()
            else: # open the image
                # if the image is a dicom image convert it to PIL Image for ease of handling
                if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                    ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                    new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                    scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                    scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                    image_array = np.asarray(scaled_image) # place the converted dicom pixel data in an array to be read like other image types

                elif self.num_channels != 1 or self.color == "Binary": # convert to greyscale if not greyscale, including binary because binary stores true and false
                    image_array = np.asarray(Image.open(self.file_path).convert("L")) # read a 2D array of image pixel data

                else:
                    image_array = np.asarray(Image.open(self.file_path)) # read a 2D array of image pixel data
            
                width = image_array.shape[1]
                height = image_array.shape[0]
                new_width = int(np.ceil(factor * width))
                new_height = int(np.ceil(factor * height))
                
                self.nearest_neighbor(image_array, height, width, new_height, new_width, factor) # call nearest neighbor interpolation function

                # call bilinear and show image
                im = Image.fromarray(np.uint8(self.resize_image(image_array, new_height, new_width, factor)))
                qimg = im.toqpixmap()
                self.bilinear_label.clear()
                self.bilinear_label.setPixmap(qimg)
                self.new_width_line.clear() 
                self.new_height_line.clear()
                self.new_width_line.setText(str(new_width)) # show new width after zooming
                self.new_height_line.setText(str(new_height)) # show new height after zooming

                # give user prompt to open the 2nd tab
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("Hello User")
                msg.setInformativeText('Please Open Tab 2 to view zoomed images')
                msg.setWindowTitle("Image Viewer")
                msg.exec_()
               
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def plot_T(self):
        # generate an image with a black background and white T with the specified dimensions
        T_image = np.zeros((128,128))
        T_image[29:49,29:99] = 255
        T_image[49:99,54:74] = 255
        # plot the array values in matplotlib to show it as an image with axes
        self.scene1.clear()    
        figure = Figure()
        ax = figure.gca() 
        ax.imshow(T_image, cmap=plt.cm.gray) 
        canvas = FigureCanvas(figure) 
        self.scene1.addWidget(canvas) 
        self.original_view.setScene(self.scene1) 

    def rotate_nearest(self):
        try:
            # 2D array to generate T image
            image = np.zeros((128,128))
            image[29:49,29:99] = 255
            image[49:99,54:74] = 255
            # get input angle from the user and multiply by negative to follow the convention that +ve angles rotate anticlockwise
            angle = float(-self.angle_spinbox.value())
            # display the direction of rotation for the input angle
            self.direction_line.clear()
            if angle > 0:
                self.direction_line.setText("Clockwise")
            elif angle < 0:
                self.direction_line.setText("Anticlockwise")
            elif angle == 0:
                self.direction_line.setText("No Rotation")
            #calculate sin and cos the angle
            angle = math.radians(angle)                               
            cosine = math.cos(angle)
            sine = math.sin(angle)
            width = image.shape[1]   
            height = image.shape[0]                               
            # array of zeros of the image size
            output = np.zeros((height,width))
            # Find the center of the image to perform fixed point rotation
            center_height = round(((height + 1)/2))    
            center_width = round(((width + 1)/2)) 

            for i in range(height):
                for j in range(width):
                    # coordinates of pixel with respect to the center of image
                    y = i - center_height                   
                    x = j - center_width 
                    
                    # coordinates of pixel with respect to the rotated image using condition for nearest neighbor interpolation (rotation matrix)
                    if x > (height - 1) or y > (width -1):
                        new_y = int(-x * sine + y * cosine)
                        new_x = int(x * cosine + y * sine)
                    else:
                        new_y = round(-x * sine + y * cosine)
                        new_x = round(x * cosine + y * sine)
                        
                    # return the center coordinates to the center of the image rather than zero
                    new_y = new_y + center_height
                    new_x = new_x + center_width
                    
                    # if the coordinate is within the box of the image place the new values into thhe output array that was filled with zeros
                    if 0 <= new_x < width and 0 <= new_y < height and new_x>=0 and new_y>=0:
                        output[i,j] = image[new_y,new_x]  

            # display the edited image
            self.scene2.clear()    
            figure = Figure()
            ax = figure.gca() 
            ax.imshow(output, cmap=plt.cm.gray) 
            canvas = FigureCanvas(figure) 
            self.scene2.addWidget(canvas) 
            self.edited_view.setScene(self.scene2) 
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def rotate_bilinear(self):
        try:
            # 2D array to generate T image
            image = np.zeros((128,128))
            image[29:49,29:99] = 255
            image[49:99,54:74] = 255
            # get input angle from the user and multiply by negative to follow the convention that +ve angles rotate anticlockwise
            angle = float(-self.angle_spinbox.value())
            # display the direction of rotation for the input angle
            self.direction_line.clear()
            if angle > 0:
                self.direction_line.setText("Clockwise")
            elif angle < 0:
                self.direction_line.setText("Anticlockwise")
            elif angle == 0:
                self.direction_line.setText("No Rotation")

            #calculate sin and cos the angle
            angle = math.radians(angle)                               
            cosine = math.cos(angle)
            sine = math.sin(angle)
            width = image.shape[1]   
            height = image.shape[0]                               
            # array of zeros of the image size
            output = np.zeros((height,width))
            # Find the center of the image to perform fixed point rotation
            center_height = round(((height + 1)/2))    
            center_width = round(((width + 1)/2)) 

            for i in range(height):
                for j in range(width):
                    #co-ordinates of pixel with respect to the center of image
                    y = i - center_height
                    x = j - center_width 
                    
                    #co-ordinate of pixel with respect to the rotated image (rotation matrix)
                    new_y = -x * sine + y * cosine
                    new_x = x * cosine + y * sine

                    # return the center coordinates to the center of the image rather than zero
                    new_y = new_y + center_height
                    new_x = new_x + center_width

                    # use bilinear interpolation to find the missing pixel values
                    output[i, j] = self.bilinear_interpolation(image, new_y, new_x)

            # display the edited image
            self.scene2.clear()    
            figure = Figure()
            ax = figure.gca() 
            ax.imshow(output, cmap=plt.cm.gray) 
            canvas = FigureCanvas(figure) 
            self.scene2.addWidget(canvas) 
            self.edited_view.setScene(self.scene2) 
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def shear(self):
        try:
            # 2D array to generate T image
            image = np.zeros((128,128))
            image[29:49,29:99] = 255
            image[49:99,54:74] = 255
            width = image.shape[1]   
            height = image.shape[0]                               
            # array of zeros same as the image size
            output = np.zeros((height,width))
            #Find the center of the image
            center_height = round(((height + 1)/2))    
            center_width = round(((width + 1)/2)) 

            for i in range(height):
                for j in range(width):
                    # co-ordinates of pixel with respect to the center of image
                    y = i - center_height               
                    x = j - center_width 
                    
                    # shearing using horizontal shearing matrix
                    new_x = x - y
                    new_y = y
                    
                    # return the center of the image to the center
                    new_y = new_y + center_height
                    new_x = new_x + center_width
                    
                    if 0 <= new_x < width and 0 <= new_y < height and new_x>=0 and new_y>=0:
                        output[i,j] = image[new_y,new_x]

            self.scene2.clear()    
            figure = Figure()
            ax = figure.gca() 
            ax.imshow(output, cmap=plt.cm.gray) 
            canvas = FigureCanvas(figure) 
            self.scene2.addWidget(canvas) 
            self.edited_view.setScene(self.scene2)
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def show_image_tab4(self):
        try:
            if self.check_tab() == 3:
                self.origin_img_lbl.clear()
                if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                    ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                    new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                    scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                    scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                    image_np_array = np.asarray(scaled_image) # place the converted dicom pixel data in an array to be read like other image types
                    self.height_img = ds.Rows
                    self.width_img = ds.Columns
                else:
                    image = Image.open(self.file_path).convert("L")
                    self.height_img = image.size[0]
                    self.width_img = image.size[1]
                    image_np_array = np.asarray(Image.open(self.file_path).convert("L"))
                    
                # display the results
                im = Image.fromarray(np.uint8(image_np_array))
                qimg = im.toqpixmap()
                if self.height_img > 600 or self.width_img > 340:
                    self.origin_img_lbl.setPixmap(qimg.scaled(600,340))
                else: 
                    self.origin_img_lbl.setPixmap(qimg)

                self.scene3.clear()
                self.scene4.clear()
                self.equal_img_lbl.clear()
        except:
            self.scene3.clear()
            self.scene4.clear()
            self.equal_img_lbl.clear()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def histogram(self, image_array):
        try:
            resolution = len(image_array)
            # count the frequencies of each pixel value
            counts = Counter(image_array) # ordered from most frequent to least
            pixel_value = np.zeros(self.max_depth)
            self.normalized_count = np.zeros(self.max_depth)
            count = np.zeros(self.max_depth) # counts is of type Counter, so convert to array to be able to divide by resolution 

            # for the values from 0 to max (255 for example)
            for i in range(self.max_depth):
                count[i] = int(counts[i]) # place values in the array
                self.normalized_count[i] = (counts[i]) / resolution # calculate the normalized count (probability), y axis of histogram
                pixel_value[i] = int(i) # array of possible pixel values as integers, x axis of histogram
            result = [pixel_value, self.normalized_count]
            return result
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def original_histogram(self):
        try:
            self.histogram_flag = 1
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                    ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                    new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                    scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                    self.max_depth = 256
                    scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                    image_np_array = np.asarray(scaled_image)
                    self.height_img = ds.Rows
                    self.width_img = ds.Columns
            else:
                # open the image and convert to greyscale
                image = Image.open(self.file_path).convert("L")
                self.height_img = image.size[0]
                self.width_img = image.size[1]
                # read image pixel values as an array
                image_np_array = np.asarray(Image.open(self.file_path).convert("L")) # read a 2D array of image pixel data
                # calculate max possible pixel values (0 - 255 for example)
                max = np.amax(image)
                depth = math.ceil(math.log((int(max) + 1), 2))
                self.max_depth = 2 ** depth

            # reshape image array to be 1D for further use of the array, note that image_1d is not writable because it's the image array 
            self.image_1d = image_np_array.reshape(-1)
        
            # relocate the values into a writable array to be used when writing values into the equalized image
            hist_input = self.histogram(self.image_1d) # call histogram function with the image array

            # display the result
            self.scene3.clear()    
            fig = Figure()
            ax = fig.gca() 
            ax.bar(hist_input[0], hist_input[1]) 
            # ax.set_ylim(top = 1)

            canvas = FigureCanvas(fig) 
            self.scene3.addWidget(canvas) 
            canvas.setGeometry(0, 0, 580, 320)
            self.orig_hist_view.setScene(self.scene3)

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
                    
    def equalized_image(self):
        try:
            if self.histogram_flag != 1: # check that buttons are pressed in the correct sequence
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("ERROR")
                msg.setInformativeText('Please click button 1 first')
                msg.setWindowTitle("Error")
                msg.exec_()
            else:
                self.histogram_flag = 2
                new_pixel_vals = np.zeros(self.max_depth) # create array for new pixel values (sk)
                CDF = 0
                for i in range(len(self.normalized_count)):
                    CDF += self.normalized_count[i] # calculate CDF for new pixel values
                    new_pixel_vals[i] = round((self.max_depth -1) * CDF) # calculate the new pixel values
                    
                equalized_image = np.zeros((self.width_img, self.height_img)) # create array of image size
                self.equalized_image_1d = equalized_image.reshape(-1) # reshape it to 1D for ease of comparison with new pixel values
                for i in range(len(self.image_1d)):
                    self.equalized_image_1d[i] = new_pixel_vals[self.image_1d[i]] # replace the new pixel value with the corresponding old pixel value
                
                # display the results
                im = Image.fromarray(np.uint8(equalized_image))
                qimg = im.toqpixmap()
                self.equal_img_lbl.clear()
                # self.equal_img_lbl.setPixmap(qimg)
                if self.height_img > 600 or self.width_img > 340:
                    self.equal_img_lbl.setPixmap(qimg.scaled(600,340))
                else: 
                    self.equal_img_lbl.setPixmap(qimg)
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def equalized_histogram(self):
        try:
            if self.histogram_flag != 2: # check that buttons are pressed in the correct sequence
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("ERROR")
                msg.setInformativeText('Please click button 2 first')
                msg.setWindowTitle("Error")
                msg.exec_()
            else:
                self.histogram_flag = 3 
                # call the histogram function on the equalized image's array
                hist = self.histogram(self.equalized_image_1d)

                # display the result
                self.scene4.clear()    
                fig = Figure()
                ax = fig.gca() 
                ax.bar(hist[0], hist[1]) 
                # ax.set_ylim(top = 1)
                canvas = FigureCanvas(fig) 
                self.scene4.addWidget(canvas) 
                canvas.setGeometry(0, 0, 580, 320)
                self.equal_hist_view.setScene(self.scene4)
        except: 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def show_image_tab5(self):
        try:
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                image_np_array = np.asarray(scaled_image)
                img_height = ds.Rows
                img_width = ds.Columns
            else:
                # open the image and convert to greyscale
                image = Image.open(self.file_path).convert("L")
                image_np_array = np.asarray(Image.open(self.file_path).convert("L"))
                img_height = image.size[0]
                img_width = image.size[1]

            im = Image.fromarray(np.uint8(image_np_array))
            qimg = im.toqpixmap()
            self.orig_img_lbl_5.clear()
            
            if img_height > 600 or img_width > 600:
                self.orig_img_lbl_5.setPixmap(qimg.scaled(600,600))
            else: 
                self.orig_img_lbl_5.setPixmap(qimg)
        except: 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def clear_median_label(self):
        self.median_img_lbl.clear()
    
    def clear_higboost_label(self):
        self.unsharp_img_lbl.clear()   
    
    def highboost_enhancement(self):
        try:
            kernel_size = int(self.kernel_size_input.value())
            multiplication_factor = float(self.mult_factor_input.value())
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                image_np_array = np.asarray(scaled_image)
                img_height = ds.Rows
                img_width = ds.Columns
            else:
                # open the image and convert to greyscale
                image = Image.open(self.file_path).convert("L")
                image_np_array = np.asarray(Image.open(self.file_path).convert("L"))
                img_height = image.size[1]
                img_width = image.size[0]

            # box kernel creation
            kernel = np.full((kernel_size,kernel_size), 1/(kernel_size ** 2)) # no need to flip convolution = correlation because kernel is symmetric
            kernel_pad = kernel_size // 2 # size to pad according to kernel size
            # image zero padding
            padded_image = zero_pad(kernel_size, img_height, img_width, image_np_array)
            padded_width = padded_image.shape[1]
            padded_height = padded_image.shape[0]

            # convolution
            convoluted_img = np.zeros((padded_height, padded_width))
            # print(convoluted_img)
            for i in range(img_height):
                for j in range(img_width):
                    sum_prod = 0 # to calculate center pixel value
                    for ii in range(kernel_size):
                        for jj in range(kernel_size):
                            sum_prod += kernel[ii][jj] * padded_image[ii + i][jj + j] # center pixel = sum of products of kernel over image
                            # we add i and to the indeces in the padded image because we won't multiply by the padding
                    convoluted_img[i + kernel_pad][j + kernel_pad] = sum_prod  # place the calculated, filtered values in a new 2D array

            subtraction_img = padded_image - convoluted_img # perform subtraction between original (padded because of size) image and filtered one
            multiplication_img  = (multiplication_factor * subtraction_img) + padded_image # perform multiplication by a factor and then add the original image to have a highboost filter
            final = rescale_intensities(multiplication_img, padded_width, padded_height) # rescale by clipping the values from 0 to 255
            clipped_image = final[kernel_pad : img_height + kernel_pad, kernel_pad : img_width + kernel_pad] # remove the border created by the padding
            # display the image in UI
            im = Image.fromarray(np.uint8(clipped_image))
            qimg = im.toqpixmap()
            self.unsharp_img_lbl.clear()
            if img_height > 600 or img_width > 600:
                self.unsharp_img_lbl.setPixmap(qimg.scaled(600,600))
            else: 
                self.unsharp_img_lbl.setPixmap(qimg)
        except: 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()
   
    def seasoning(self):
        try:
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                image_np_array = np.asarray(scaled_image)
                img_height = ds.Rows
                img_width = ds.Columns
            else:
                # open the image and convert to greyscale
                image = Image.open(self.file_path).convert("L")
                image_np_array = np.asarray(Image.open(self.file_path).convert("L"))
                img_height = image.size[1]
                img_width = image.size[0]

            resolution = img_width * img_height
            percentage = int(self.percent_noise_input.value()) # percentage of pixels the user wants to add noise to

            self.seasoned_img = np.zeros((img_height, img_width)) # self.seasoned_img is a writable version of image_np_array
            for i in range(img_height):
                for j in range(img_width):
                    self.seasoned_img[i][j] = image_np_array[i][j]

            # seasoning 
            number_of_pixels = int((percentage/200) * resolution) # number of pixels to be changed in EACH color (black and white)

            for i in range(number_of_pixels):
                # randomly choose pixels to change their value
                y_coord=random.randint(0, img_height - 1) 
                x_coord=random.randint(0, img_width - 1)
                self.seasoned_img[y_coord][x_coord] = 255 # change their values to white

            for i in range(number_of_pixels):
                # randomly choose pixels to change their value
                y_coord=random.randint(0, img_height - 1)
                x_coord=random.randint(0, img_width - 1)
                self.seasoned_img[y_coord][x_coord] = 0 # change their values to black

            im = Image.fromarray(np.uint8(self.seasoned_img))
            qimg = im.toqpixmap()
            self.salt_pep_img_lbl.clear()
            if img_height > 600 or img_width > 600:
                self.salt_pep_img_lbl.setPixmap(qimg.scaled(600,600))
            else: 
                self.salt_pep_img_lbl.setPixmap(qimg)
        except: 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def median_filter(self):
        try:
            kernel_size = int(self.kernel_size_in_snp.value())
            kernel_pad = kernel_size // 2 # size to pad according to kernel size
            img_width = self.seasoned_img.shape[1]
            img_height = self.seasoned_img.shape[0]
            # image zero padding
            padded_image = zero_pad(kernel_size, img_height, img_width, self.seasoned_img)
            padded_width = padded_image.shape[1]
            padded_height = padded_image.shape[0]
            # set the array for filtered image and kernel size (user input)
            median_img = np.zeros((padded_height, padded_width))
            window = np.zeros((kernel_size, kernel_size))

            for i in range(img_height):
                for j in range(img_width):
                    for ii in range(kernel_size):
                        for jj in range(kernel_size):
                            window[ii][jj] = padded_image[ii + i][jj + j] # get current window
                            widow1d = window.flatten() # flatten to handle later
                            sorted_arr = quicksort(widow1d) # sort to find median
                            median_val = sorted_arr[int(np.floor(((kernel_size **2) /2)))] # get median term
                
                    median_img[i + kernel_pad][j + kernel_pad] = median_val  # store in new image

            clipped_image = median_img[kernel_pad : img_height  + kernel_pad, kernel_pad : img_width + kernel_pad] # remove padding
            im = Image.fromarray(np.uint8(clipped_image))
            qimg = im.toqpixmap()
            
            if img_height > 600 or img_width > 600:
                self.median_img_lbl.setPixmap(qimg.scaled(600,600))
            else: 
                self.median_img_lbl.setPixmap(qimg)
        except: 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def show_image_tab7(self):
        try:
            self.original_img_view.clear()
            if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                image_np_array = np.asarray(scaled_image) # place the converted dicom pixel data in an array to be read like other image types
                self.height_img = ds.Rows
                self.width_img = ds.Columns
            else:
                image = Image.open(self.file_path).convert("L")
                self.height_img = image.size[0]
                self.width_img = image.size[1]
                image_np_array = np.asarray(Image.open(self.file_path).convert("L"))
                
            # display the results
            im = Image.fromarray(np.uint8(image_np_array))
            qimg = im.toqpixmap()
            if self.height_img > 600 or self.width_img > 600:
                self.original_img_view.setPixmap(qimg.scaled(600,600))
            else: 
                self.original_img_view.setPixmap(qimg)

        except:

            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('An error has occured')
            msg.setWindowTitle("Error")
            msg.exec_()

    def fourier_transform(self):
        if magic.from_file(self.file_path) == 'DICOM medical imaging data' or magic.from_file(self.file_path) == 'TIFF image data, little-endian': # call dicom function if file is dicom by checking file type first
                    ds = pydicom.dcmread(self.file_path) # read the dicom file dataset
                    new_image = ds.pixel_array.astype(float) # read the pixel array values as floats
                    scaled_image = (np.maximum(new_image, 0) / new_image.max()) * 255.0 # scale the values in 8 bits with values between 0 and 255
                    scaled_image = np.uint8(scaled_image) # 8-bit unsigned integer, used for arrays representing images with the 3 color channels having small integer values (0 to 255).
                    image_np_array = np.asarray(scaled_image) # place the converted dicom pixel data in an array to be read like other image types

        else:
            image_np_array = np.asarray(Image.open(self.file_path).convert("L"))

        ft_image_np = np.fft.fft2(image_np_array)
        fshift = np.fft.fftshift(ft_image_np)
        real_component = fshift.real
        imaginary_component = fshift.imag

        self.magnitude = np.sqrt((real_component ** 2) + (imaginary_component ** 2))
        self.phase = np.arctan2(imaginary_component, real_component)
        self.log_magnitude = np.log(self.magnitude + 1)
        self.log_phase = np.log(self.phase + 2 * math.pi)

        self.display_magnitude()
        self.display_log_magnitude()
        self.display_phase()
        self.display_log_phase()

    def display_magnitude(self):

        self.scenem.clear()
        # show the dicom image in the GUI       
        figure = Figure() # using matplotlib's Figure which holds all plot elements
        ax = figure.gca() # get the current axes
        # remove appearance of axes to show the image
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.imshow(self.magnitude, interpolation = "None", cmap="gray", vmin = 0) # take pixel data from dicom file and display the data as an image in the bone color map
        canvas = FigureCanvas(figure) # pass the figure to Figure canvas which handles all the details of talking to user interface toolkits (pyqt in this case)
        self.scenem.addWidget(canvas) # add a widget to the scene in QGraphicsView and use canvas to interface with matplotlib
        self.magnitude_view.setScene(self.scenem) # set the scene in the ui's graphics view
        self.magnitude_view.fitInView(self.scenem.sceneRect())

    def display_log_magnitude(self):

        self.sceneml.clear()
        # show the dicom image in the GUI       
        figure = Figure() # using matplotlib's Figure which holds all plot elements
        ax = figure.gca() # get the current axes
        # remove appearance of axes to show the image
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.imshow(self.log_magnitude, interpolation = "None", cmap="gray", vmin = 0) # take pixel data from dicom file and display the data as an image in the bone color map
        canvas = FigureCanvas(figure) # pass the figure to Figure canvas which handles all the details of talking to user interface toolkits (pyqt in this case)
        self.sceneml.addWidget(canvas) # add a widget to the scene in QGraphicsView and use canvas to interface with matplotlib
        self.magnitude_log_view.setScene(self.sceneml) # set the scene in the ui's graphics view
        self.magnitude_log_view.fitInView(self.sceneml.sceneRect())

    def display_phase(self):

        self.sceneph.clear()
        # show the dicom image in the GUI       
        figure = Figure() # using matplotlib's Figure which holds all plot elements
        ax = figure.gca() # get the current axes
        # remove appearance of axes to show the image
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.imshow(self.phase, interpolation = "None", cmap="gray", vmin = 0) # take pixel data from dicom file and display the data as an image in the bone color map
        canvas = FigureCanvas(figure) # pass the figure to Figure canvas which handles all the details of talking to user interface toolkits (pyqt in this case)
        self.sceneph.addWidget(canvas) # add a widget to the scene in QGraphicsView and use canvas to interface with matplotlib
        self.phase_view.setScene(self.sceneph) # set the scene in the ui's graphics view
        self.phase_view.fitInView(self.sceneph.sceneRect())
        
    def display_log_phase(self):

        self.scenephl.clear()
        # show the dicom image in the GUI       
        figure = Figure() # using matplotlib's Figure which holds all plot elements
        ax = figure.gca() # get the current axes
        # remove appearance of axes to show the image
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.imshow(self.log_phase, interpolation = "None", cmap="gray", vmin = 0) # take pixel data from dicom file and display the data as an image in the bone color map
        canvas = FigureCanvas(figure) # pass the figure to Figure canvas which handles all the details of talking to user interface toolkits (pyqt in this case)
        self.scenephl.addWidget(canvas) # add a widget to the scene in QGraphicsView and use canvas to interface with matplotlib
        self.phase_log_view.setScene(self.scenephl) # set the scene in the ui's graphics view
        self.phase_log_view.fitInView(self.scenephl.sceneRect())

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()
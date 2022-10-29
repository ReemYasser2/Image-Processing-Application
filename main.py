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
  

class Ui(QtWidgets.QMainWindow): # class to load .ui in constructor
    def __init__(self):
        super(Ui, self).__init__() # call the inherited classes __init__ method
        uic.loadUi('image_viewer.ui', self) # load the .ui file
        self.setWindowTitle("Image Viewer")
        self.show() # Show the GUI
        self.scene = QGraphicsScene()
        self.image_view.setScene(self.scene) #set a graphics scene in QGraphicsView
        self.scene1 = QGraphicsScene()
        self.original_view.setScene(self.scene1) #set a graphics scene in QGraphicsView
        self.action_open.triggered.connect(self.browse_files) # call browse_files when Open is clicked in menubar
        self.zoom_button.clicked.connect(self.zooming)
        self.plot_T()
    
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
            self.width_line.clear() 
            self.height_line.clear() 
            self.size_line.clear()
            self.color_line.clear()
            self.depth_line.clear()
            self.modality_line.clear()
            self.name_line.clear()
            self.age_line.clear()
            self.part_line.clear()
            self.scene.clear()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("ERROR")
            msg.setInformativeText('Please choose a file to open')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def show_image(self): # shows any common image type (jpg, bmp, png, etc.)
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
            self.width_line.clear() 
            self.height_line.clear() 
            self.size_line.clear()
            self.color_line.clear()
            self.depth_line.clear()
            self.modality_line.clear()
            self.name_line.clear()
            self.age_line.clear()
            self.part_line.clear()

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
            print(ds.pixel_array.shape)
        
            # clearing past displayed data
            self.width_line.clear() 
            self.height_line.clear() 
            self.size_line.clear()
            self.color_line.clear()
            self.depth_line.clear()
            self.modality_line.clear()
            self.name_line.clear()
            self.age_line.clear()
            self.part_line.clear()

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

    def nearest_neighbor(self, image_array, new_height, new_width, factor):
        new_arr = np.zeros((new_height, new_width)) # 2D array of zeros with size same as zoomed image
        # replace the zeros with the correct pixel data 
        for i in range(new_height):
            for j in range(new_width):
                new_arr[i,j] =  image_array[int(np.floor(i/factor)), int(np.floor(j/factor))] 

        # display image
        im = Image.fromarray(np.uint8(new_arr))
        qimg = im.toqpixmap() # create qpixmap to display image in Qlabel
        self.neighbor_label.clear() # clear past images
        self.neighbor_label.setPixmap(qimg) # set the created pixmap in the label
        # self.neighbor_label.setPixmap(qimg.scaled(650,700)) # for re-scaling
    
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
                new_width = int(factor * width)
                new_height = int(factor * height)
                
                self.nearest_neighbor(image_array, new_height, new_width, factor) # call nearest neighbor interpolation function

                # call bilinear and show image
                im = Image.fromarray(np.uint8(self.resize_image(image_array, new_height, new_width, factor)))
                qimg = im.toqpixmap()
                self.bilinear_label.clear()
                self.bilinear_label.setPixmap(qimg)
                # self.bilinear_label.setPixmap(qimg.scaled(650,700)) # for re-scaling
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
        T_image = np.zeros((128,128))
        T_image[29:49,29:99] = 255
        T_image[49:99,54:74] = 255
        self.scene1.clear()    
        figure = Figure()
        ax = figure.gca() 
        ax.imshow(T_image, cmap=plt.cm.gray) 
        canvas = FigureCanvas(figure) 
        self.scene1.addWidget(canvas) 
        self.original_view.setScene(self.scene1) 

    

            
app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()


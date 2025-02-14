# Optical-sectioning confocal structured illumination microscopy (OS-CSIM)

# Overview

We introduce the confocal imaging principle into the structured illumination microscopy and propose ann optical-sectioning confocal structured illumination microscopy to improve the signal-to-noise ratio and the imaging depth of fluorescent optical sectioning imaging [[paper]](https://opg.optica.org/oe/fulltext.cfm?uri=oe-32-18-32550&id=554913). To realize confocal imaging in structured illumination microscopy, we firstly calibrate the conjugate relationship between the spatial light modulator and the camera, then use the conjugate relationship to generate digital pinhole. Next, we use each pixel of the camera as a single-pixel detector and use the principles of dual imaging and Fourier single-pixel imaging to reconstruct dual image to separate the  signal from the conjugate object point and the noise from the non-conjugated object points. Third, we use the generated digital pinhole to extract the conjugated signal from the dual imagge and put it back to the camera pixel. Some Fourier single-pixel reconstruction programs refer to the open-source code provided by Zhang Zibang [[code]](https://github.com/zibangzhang/Fourier-single-pixel-imaging). The whole program mainly consists of three modules:

(1) Generate projection pattern; 

(2) Calibrate the congjuate relationship between the spatial light modulator and the camera; 

(3) Reconstruct optical-sectioning image. 

In this study, we also made an important discovery: there’s no need to reconstruct the dual image. Only by relying on the measured Fourier coefficients and the digital pinhole can we calculate the conjugate signal. Therefore, we have provided two reconstruction programs. One of them requires the use of dual images when reconstructing light slice images, while the other does not require the use of dual images during the reconstruction of light slice images. Using the latter reconstruction program (i.e., the one that does not require the use of dual images) can significantly reduce the time required for image reconstruction and effectively improve the reconstruction efficiency.

We have also provided a program for reconstruction optical sectioning image with OS-SIM.

# raw data in this study can be downloaded with the following link: 
[https://drive.google.com/drive/folders/1zCXlDlkdB2lWyqny-jmcRhDHQwoql1FR?usp=sharing](https://drive.google.com/drive/folders/15WXBXKNKnL3vGFGChvT-2vV-0_vHCidb?usp=sharing)



# Citation

If you use this code and relevant data, please cite the corresponding paper where original methods appeared:

Weishuai Zhou, Manhong Yao, Xi Lin, Quan Yu, Junzheng Peng, and Jingang Zhong, "Confocal structured illumination microscopy for improving the signal-to-noise ratio and depth of fluorescent optical section imaging," Opt. Express 32, 32550-32563 (2024). https://doi.org/10.1364/OE.536711.

# Correspondence

Should you have any questions regarding this code and the corresponding results, please contact Junzheng Peng (junzpeng@jun.edu.cn). 




Team 2a of the Big Data REU program at UMBC for Summer 2026

**Title:** Experimental Prompt Gamma Scatter Classification and Reorientation for Proton Beam Image Reconstruction

**Team Members:** Ayan Kabaria<sup>1</sup>, Muhammad Khalid<sup>2</sup>, Sophia Lopez<sup>3</sup>, Abby Nam<sup>4</sup>, Ertan Dogan<sup>5</sup>, Sidhya Pathak<sup>6</sup>, Victor Sandrin<sup>7</sup>, Xueying Sun<sup>8</sup>, Ehsan Shakeri<sup>9</sup>, Harrison Lewis<sup>9</sup>, Hussam Fateen<sup>9</sup>, Matthias K. Gobbert<sup>9</sup>, Farshad Safavi<sup>10</sup>, Ananta Chalise<sup>10</sup>, Lei Ren<sup>11</sup>, Stephen W. Peterson<sup>12</sup>, and Jerimy C. Polf<sup>13</sup>

<sup>1</sup> River Hill High School, Howard County, Maryland  
<sup>2</sup> Department of Mathematics, Baruch College, City University of New York  
<sup>3</sup> Department of Mathematics and Statistics, University of North Carolina at Greensboro  
<sup>4</sup> Department of Psychology, Lafayette College  
<sup>5</sup> A. James Clark School of Engineering, University of Maryland, College Park  
<sup>6</sup> Department of Computer Science, University of Virginia  
<sup>7</sup> Department of Neuroscience, University of Arizona  
<sup>8</sup> Department of Information Technology and Management, Illinois Institute of Technology  
<sup>9</sup> Department of Mathematics and Statistics, University of Maryland, Baltimore County  
<sup>10</sup> Department of Radiation Oncology, University of Maryland School of Medicine  
<sup>11</sup> Department of Radiation Oncology, Northwestern University  
<sup>12</sup> Department of Physics, University of Cape Town, South Africa  
<sup>13</sup> M3D, Inc.

**Abstract:** Prompt gamma rays emitted by cells during proton beam radiotherapy can be detected by Compton cameras and used to reconstruct images to verify accurate dose delivery and Bragg peak location. However, image reconstruction is complicated by misclassified scattering events because of the limited temporal resolution of the cameras. To address this, machine learning models were trained on simulated datasets generated using Geant4 and MCDE and implemented through the Big Data REU Integrated Development and Experimentation (BRIDE) platform to classify and filter scatter events for more accurate image reconstruction. After identifying the best-performing models, predictions using a deep Fully Connected Network (FCN) were created for 24 cases of experimental phantom triples data. Images were reconstructed on CORE using Kernel Expectation-Maximization (KEM) and Simple Back Projection (SBP) algorithms on both raw and repaired data, then compared to assess the accuracy of the Bragg peak estimation and noise reduction. SBP reconstructions had more complete beam reconstructions and Bragg peak profiles than KEM reconstructions, however still included noise in 2-D beam profiles. KEM 1D profile reconstructions improved after a reorientation script was implemented to transform the scatter data before classification, which reoriented the data such that it was consistent with training, and reversed before reconstruction.

### Navigating the Repository
The names of each directory in this repository correspond to a step in the data generation pipeline. *1-PJMC* and *2-MCDE* are integral steps for simulating scatter data and applying detector effects to make the data more realistic, respectively, but were not used for the results of our study. *3-BRIDE* contains the files and scripts used to apply machine learning models to classify the scatters as well as specific instructions for using BRIDE for experimental triples data. *4-CORE* creates 1-dimensional and 2-dimensional image reconstructions of the data, and also includes specific instructions for applying the scripts to experimental triples data. 

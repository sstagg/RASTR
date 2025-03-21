# RASTR
Repository for reconstruction of average subtracted tubular regions(RASTR) code.
Also contains code for  SPOT-RASTR, diameTR and helpful python scripts for cryo-EM data processing.

## Installation
### Requirements
- python >= 3.6
- cupy, numpy, scipy, matplotlib, tkinter, pandas, mrcfile, setuptools
- Relion >= 3.0
- We recommend create a conda environment for RASTR
    ```bash
    conda create --name RASTR python=3.9 scipy cupy numpy matplotlib tk pandas mrcfile setuptools -c conda-forge
    ```

### installation
- Clone the code
    ``` bash
    git clone https://github.com/sstagg/RASTR.git
    
    cd RASTR

    conda activate RASTR

    pip install -e .
    ```



## Usages


<details>
<summary>RASTR</summary>
A star file containing particle information is used as input, represented by 'particles.star'. During the process, two parameter optimization windows will pop up. If the default values yield poor performance, refer to 'diameTR' for optimization guidance.

- 1 To get the psi angles first,
    ```bash
    # Go the main path
    cd /path/for/particles/star
    diameTR --i particles.star --o particles_p -p
    ```
- 2 Determine shift and diameter.
    ```bash
    diameTR --i particles_p.star --o particles_pds -d -s
    ```
    Diameter distribution window will pop up for you to thresholding the diamters.

- 3 Create azimuthal average model
    ```bash
    # Assign values for tilt, rot angles
    changestar --i particles_pds.star --o particles_pds.star -rot r360 -tilt 90

    # Reconstruct
    relion_reconstruct --i particles_pds.star --o particles.mrc --ctf

    # Average along y
    azavg particles.mrc
    ```

- 4 Get correct weighing

    In Relion GUI, choose 3D classification, change below parameters.

    Input images STAR file: particles_pds.star;
    Reference 'particlesazavg.mrc';
    Number of classes: 1;
    Number of iterations: 5; 
    Perform image alignment: No;

    Run!

- 5 Create mask
    ```bash
    # Create a sphere mask
    createmask  boxsize center radius pixel_size
    ```

- 6 Create RASTR particles
    ```bash
    # Go to the Relion 3Dclass job path
    azavg run_it005_class001.mrc
    
    # Go back to main path
    cd ../..

    # Run RASTR. Replace pixel_size, rootname, spheremask.mrc with correct filename
    RASTR --star_in Class3D/job001/run_it005_data.star  --model run_it005_class001azavg.mrc  --angpix pixel_size  -k -o rootname -ma spheremask.mrc  -al 0,90,180,270  --pad 3
    ```
<br>
</details>



<details>
<summary>SPOT-RASTR</summary>
<br>
- Have your particle stack and star file ready.

- Follow the same steps of 1-6 of RASTR to create subtracted particles

- Remove duplicates from overlapping filament particle picking
    ```bash
    Under construction.
    ```

</details>






<details>
<summary>diameTR</summary>
<br>
Have your particle stack and star file ready.

- Determine psi angles
    ```bash
    diameTR --i particles.tar --o particles_p -p
    ```

    A optimiser window will pop up.
    Press R to navigator random slice. Monitor the particle image and the bottom 1D curves.
    A good set of parameters should have filament particle horizontally oriented and the 1D curves with two clear peaks.
    In case default parameters fails, perform following optimization.

    1. Adjust length. Find the high intensity line crossing the center. Estimate the length. Put a slightly bigger number for length.
    For FT strategy the optimal length is usuall small. For AC strategy, the lenght is big.
    2. Sigma (int). This is the width of Gaussian filter. The bigger sigma, the greater of low filter. Sequentially increase it and check performance.
    3. Pad. This only affect FT strategy. Set a bigger pad will increase sample rate in fourier space and thus increasing accuracy. Big pad will slower down the compution. 
    4. Other parameters don't affect accuracy.
    5. Navigate at least 50 particles to confirm the accuracy.
    

- Determine diameter
    ```bash
    diameTR --i particles_p.star --o particles_pd -d
    ```
    A optimiser window will pop up. Press R to navigator random slices. Monitor the 1D projection curve and two scatter points. A good set of parameters should have scatter points at the edge of tubules. 
    1. Change sigma. Start from 2 and sequentially increase. Usually sigma between 3 to 5 will work.
    2. Min gap. 0 is good in most cases. When sometimes the biggest peak difference appear in the middle of tubules. Based on the diameter reported in logs, put a number around 0.7 * diameter for min_gap to eliminate errors.
    3. Other parameters shall remain unchanged.
    4. Navigate at least 50 particles to confirm the accuracy.

    A window will pop up when completed for you to thresholding the diameters.

- Other options.

    -s --shift. Whether or not to compute shifts. Default False. Useful for RASTR and SPOT-RASTR.

    --classify  Thresholding diameters again. Useful when you select a big group of particles first, then separate them into smaller groups

    --particle_number Select a random small subset of particle to test performace.

    --showaverage Average all particles together with psi rotation and centering. Usefull for RASTR and SPOT-RASTR to examine accuracy.

    --average_power_spectrum Average power spectrums for helical indexing.

</details>

<details>
<summary>Other tools</summary>
<br>
csexport.py --- A wrapper of csparc2star.py in pyem. Create softlinks to make particle stack suffix as mrcs and path handling.

- Navigator to the exported job directory. 

    ```bash
    csexport J2_particles_exported.cs J2_particles_exported.star
    ```

averagefft --- Standalone script to compute averaged power spectrum of tubular images for helical indexing.

changestar --- Star file handler. Used in house to manipulate orientations, shifts, and substitute with other star files.

</details>

## References
RASTR
- Randolph, P. S., & Stagg, S. M. (2020). Reconstruction of Average Subtracted Tubular Regions (RASTR) enables structure determination of tubular filaments by cryo-EM. Journal of Structural Biology: X, 4(February), 100023. https://doi.org/10.1016/j.yjsbx.2020.100023

SPOT-RASTR
- Esfahani, B. G., Randolph, P. S., Peng, R., Grant, T., Stroupe, M. E., & Stagg, S. M. (2024). SPOT-RASTR—A cryo-EM specimen preparation technique that overcomes problems with preferred orientation and the air/water interface. PNAS Nexus, 3(8), 284. https://doi.org/10.1093/pnasnexus/pgae284

diameTR
- Peng, R., Elkhaligy, H., Grant, T., Stagg, S.M. (2025). DiameTR: A cryo-EM tool for diameter sorting of tubular samples. Submitted.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Contact
For any questions or feedback, please contact Scott M. Stagg: sstagg@fsu.edu
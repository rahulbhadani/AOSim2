AOSim2: MATLAB-based paraxial wave optics, random media, and Adaptive Optics (AO) simulation.

Most of this code is the property of Johanan L. Codona, Ph.D.

I have included some code from the MathWorks web site to compute Zernike polynomials and Mark Milton's code to compute disk harmonics.  
I am trying to clean up this repo wrt licenses.  If in doubt, ask.  
See the license file for more details.

To install, download or clone the AOSim2 repository.  It will make a subdirectory called AOSim2 with another subdir called AOSim2 (the first of many IQ tests).

Then, in MATLAB, add the following paths:
{unpack_dir}/AOSim2/AOSim2	    % The AOSim2 classes.  <br>
{unpack_dir}/AOSim2/AOSim2/utils    % The AOSim2 utility functions.
{unpack_dir}/AOSim2/AOSim2/demos    % Optional demos, mostly examples.
{unpack_dir}/AOSim2/AOSim2/data     % Don't include in path, just know that you might need things in here for demos.
{unpack_dir}/AOSim2/contrib	    % Alex Rodack's tutorial files and things made by users.  (Hopefully more in the future.)
{unpack_dir}/AOSim2/examples	    % Contributed AOSim2 examples.  I will merge my demos files in here eventually.

You really only need the AOSim2 subdir and AOSim2/utils to make this go.  The rest is for training.

I will eventually provide a more complete manual and tutorial.  I don't include a GUI because I am a strong believer in 
the power of a consise descriptive language over gratuitous pictograms. (http://catb.org/esr/writings/unix-koans/gui-programmer.html)

I am developing this as my own personal tool as my needs evolve.  I have used it to model many interesting physical situations that would 
have been very difficult otherwise.  I have also used it to explore physical situations that have led to important theoretical insights.

I also recently have looked at what I can do with AOSim2 in comparison with commercial software from dedicated wave optics modeling 
companies, and I think I can do just about anything they can.  Even more in some cases.  The best part is due to my software design that 
makes MATLAB code look like operator mathematics.  This sometimes triggers warnings from MATLAB, but makes the code very concise and readable.
A complete AO simulation with an extended turbulent atmosphere is only a few lines of code.  I wrote it as a tool for exploring physics, and 
that has been a major success.  

Unlike companies that sell software with similar capabilities, I am making this open source in the hope that it will be useful to others.
If you want to contribute code, it would be welcomed.  I do insist that the spirit of the design be followed and that objects only "know"
what they would have access to in a physical context.  For example, an AOField have a wavelength, but not an AOAperture, unless it referred 
to a stop/pass-band.  Code written from a physics perspective should be fine.

I am adding things as I work on projects.  I recently added support for GPUs, in the most non-invasive way.  As always, it helps to know 
what you're doing, but if you create an AOField F, just call F.useGPU(true) to get started.  All AOGrids and their children have the support.

I think I could pretty much make a living off using this software at this point, which is remarkable.  I would be happy to see it do the same for 
others, use the features of github to start discussions or ask questions.

I am currently working in the branch jlc1.4.

Johanan (John) Codona
Tucson, AZ.

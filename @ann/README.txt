        
        MATLAB CLASS FOR COMPUTING APPROXIMATE NEAREST NEIGHBORS
        
        WRAPPER FOR DAVID M. MOUNT & SUNIL ARYA ANN LIB


   This wrapper for Matlab was written by Shai Bagon (shai.bagon@weizmann.ac.il).
   Department of Computer Science and Applied Mathmatics
   Wiezmann Institute of Science
   http://www.wisdom.weizmann.ac.il/~bagon

        The core cpp application was written by David M. Mount and Sunil Arya
	(available from http://www.cs.umd.edu/~mount/ANN/):

	ANN is a library for approximate nearest neighbor searching,
	based on the use of standard and priority search in kd-trees
	and balanced box-decomposition (bbd) trees. Here are some
	references to the main algorithmic techniques used here:

		kd-trees:
			Friedman, Bentley, and Finkel, ``An algorithm for finding
				best matches in logarithmic expected time,'' ACM
				Transactions on Mathematical Software, 3(3):209-226, 1977.

		Priority search in kd-trees:
			Arya and Mount, ``Algorithms for fast vector quantization,''
				Proc. of DCC '93: Data Compression Conference, eds. J. A.
				Storer and M. Cohn, IEEE Press, 1993, 381-390.

		Approximate nearest neighbor search and bbd-trees:
			Arya, Mount, Netanyahu, Silverman, and Wu, ``An optimal
				algorithm for approximate nearest neighbor searching,''
				5th Ann. ACM-SIAM Symposium on Discrete Algorithms,
				1994, 573-582.

==============================================================================================

 INSTALLATION

 1. Extract all files from archive into a directory (i.e., ~/matlab_tools)
 2. Make sure you have thefollowing folders created
    ~/matlab_tools/ann_wrapper
    ~/matlab_tools/ann_wrapper/@ann
    ~/matlab_tools/ann_wrapper/@ann/private
 3. Open Matlab
 4. Add ~/matlab_tools/ann_wrapper to your path 
 5. If you haven't use Matlab's compiler yet, you need to set it up:
    >> mex -setup
    Matlab detects any compilers installed on your machine and let you choose from them.
    I suggest using visual studio on Windows machines and gcc on Linux.
 6. Go to the ann_wrapper directory
    >> cd ~/matlab_tools/ann_wrapper
 7. Compile the cpp library
    >> ann_class_compile
 8. Test that everythig is properly installed
    >> test_ann_class
 9. If there are no error messages - you are good to go.
 10. For help on the class type:
    >> doc ann

==============================================================================================

 USING THE WRAPPER

 ANN class

 Usage:
   anno = ann(pts)

 Input:
   pts - (d)x(N) matrix of d-dimensional veecotrs representing N points

 Output:
   anno - object handle to be used in class' methods

 Methods:
   ksearch:
       k-nn search
       [idx dst] = ksearch(anno, pt, k, eps, [asm])

       idx - indices of matching points:  (k)x(#points)
       dst - square distance of the points

       anno - active ann object handle
       pt   - query point(s): a (d)x(#points) matrix
       k    - required k aprox. nearest neighbors
       eps  - accuracy multiplicative factor
       asm  - allow self match flag (optional, default: true)

   prisearch:
       priority search
       [idx dst] = ksearch(anno, pt, k, eps, [asm])

       see ksearch for inputs/outputs explaination

   frsearch:
       fixed radius search
       [idx dst inr] = frsearch(anno, pt, r, k, eps, [asm])

       idx - indices of matching points
       dst - square distance of the points
       inr - how many points are in r (range query)

       anno - active ann object handle
       pt   - query point (only one this time, (d)x1 vector)
       r    - search radius (_not_ squared)
       k    - required k aprox. nearest neighbors. Will return at most k nieghbors
              of distance at most r.
       eps  - accuracy multiplicative factor
       asm  - allow self match flag (optional, default: true)

       note that numel(idx) and numel(dst) = min(k, inr)
       inr might be larger than k reflecting how many points are actually
       in the search radius.
       
   close:
      you MUST explicitly close the ann class handle when you are done with it.
      This is crucial for releasing memory and other resources.
      
      anno = close(anno);

==============================================================================================

   This software is provided under the provisions of the
   Lesser GNU Public License (LGPL).
   see: http://www.gnu.org/copyleft/lesser.html.
   This software can be used only for research purposes, you should cite
   the aforementioned papers in any resulting publication.

   The Software is provided "as is", without warranty of any kind.

==============================================================================================

   CITATION
   
   If you use this code in a scientific project, you should cite the following works
   in any resulting publication:

@INPROCEEDINGS{Mount1993,
  author = {Sunil Arya and David M. Mount},
  title = {Approximate Nearest Neighbor Queries in Fixed Dimensions},
  booktitle = {Proc. 4th Ann. ACM-SIAM Symposium on Discrete Algorithms (SODA)},
  year = {1993},
  pages = {271-280},
  owner = {bagon},
  timestamp = {2009.02.04}
}

@ELECTRONIC{Mount2006,
  author = {David M. Mount and Sunil Arya},
  month = {August},
  year = {2006},
  title = {ANN: A Library for Approximate Nearest Neighbor Searching},
  note = {version 1.1.1},
  url = {http://www.cs.umd.edu/~mount/ANN/},
  owner = {bagon},
  timestamp = {2009.02.04}
}

@ELECTRONIC{Bagon2009,
  author = {Shai Bagon},
  month = {February},
  year = {2009},
  title = {Matlab class for ANN},
  url = {http://www.wisdom.weizmann.ac.il/~bagon/matlab.html},
  owner = {bagon},
  timestamp = {2009.02.04}
}

==============================================================================================




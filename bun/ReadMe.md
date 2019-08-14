The program draws a cat on a rice bowl. The image is made up of 10 bezier patches. Two for the face, two for the upper part of the bowl, two for the stand of the bowl and four for the ears. The program also uses three different texture maps which are the following: One for the bowl, one for the face of the cat and one for the rest of the cat’s body. The body of the cat-bun has different material properties than the bowl to stress the shininess of the bowl versus the shininess of the cat. 

Copy and paste backBun.bmp, bun.bmp and cup.bmp from the bun directory into the Debug directory to run bun.exe as a stand-alone executable.

To control the mesh of the entire rendering:
*   Press "U" to increase the number samples in U direction.
*   Press "u" to decrease the number samples in U direction.
*   Press "V" to increase the number samples in V direction.
*   Press "v" to decrease the number samples in V direction.
To control the animation:
*   Press the "a" key to toggle the animation off and on.
*   Press the "s" key to perform a single step of the animation.
*   The left and right arrow keys controls the	rate of rotation around the y-axis.
*   The up and down arrow keys increase and decrease the rate of		rotation around the x-axis. 
*   In order to reverse rotational direction you must zero or reset the patch ("0" or "r").
*   Press the "r" key to reset back to initial	position, with no rotation.
*   Press "0" (zero) to zero the rotation rates.

 Here’s a picture of the program: 
 
![Front](https://dl.boxcloud.com/api/2.0/internal_files/507675752282/versions/537557209082/representations/png_paged_2048x2048/content/1.png?access_token=1!luXHLt3ee7xmPMjo8OkIvRDHklCYB5-SWC-6svKSp7YxUb7KtupH9GddnzXVOJgNNRHgecG5ihGMNTRCfmGwrbeUVs9Iq7xYEl5uCcO_mvoRHXsyloVX4H2dljJ9AGSm-iNPVWRs-6OmHhDvp8cILES-QX4K5eVlQKkVmOTTv5w1YWBcqBVQJk_6Kul73J2-CAUBmNS7I8el6GLOrvEdTHrreojja3cD1EqEqHxwS7BSui4qOHioPfVQ8mOQUU7zWrS4ggcJT4VvkXSBp_uVkJmmkQLy1HBsUFKBSIfOWZDpuB3-opokjJdmtj0BGiL5P5T-LRfRVDseurhVG9tVDMw1F0ie9_4oWZGkYlLglwlsGQafpCvToHCrCqgGFYJIFp7-mdgCZYi16Gs77eoRwrVvnMAoPytk08I2KYzD9a_DAFJ2Mf6s_8H5kzrTIxWrHbS-auFqrcjhgRd5fz7rGFKJS6Lky3FkS3isqYw7QzfYyj-EQi4GRdDjIWAQsDib05UhMwr7GxPSe9tOCwlOITkG3_XosfWOg2jNd40z9SbEDy31ms_kY4O68_wdmuMQiQ..&box_client_name=box-content-preview&box_client_version=2.14.1)

![Back of project](https://dl.boxcloud.com/api/2.0/internal_files/507663953384/versions/537545382584/representations/png_paged_2048x2048/content/1.png?access_token=1!wE8sPUtrA4U4yriZrRiUSg2H-YjEnKGUMNDFGoO4Dmt8nlldhL-K73vQxrEGKFJNG11S_aidm7AKNWJc7NK_BRgVyejLEK10RH0_a-QcUuJTP86oe6XqQ3PMF2Nu9yQsEzN_dt2sQM8IyR4rAYqPglRY9FSkg7hlQXbT9P-Th2OKh6FtFFS3w_4e7CGu5Y19u-ANg2KDJE5EVRApeYi_c1DGeuEi87QEpnJAkGQj8z15vT9XOFDSpLRipkbc5k1lhf4GZ-2T_aQ7S3CZqfD1VVNfRngRXN6xcrNGsJDAwaX6HNY8tQtAWLJt221QAcNF9RQNiBQ0JQmYAWyzra1r3A9yZF8F8sMOrJdlvet43VfJaSnLYtiCx5G1ZqcW8jk7tCOluUmth94zjNJGSrjXHDFjjc5MIu8nvfF-a6IKZZPZhngTmhMDama_ehxD6V9isrWC8aondruNotQBgi9eZ1kCZxRESK0_mKy5ty5xkMD51iZIUvEkV4yuBvoEVq6oG41LDmjlTXfbORZwx30hlWRpK35Td8Ic0eEovj61WjLd7FHGa7Hcpn673t7op_GwIg..&box_client_name=box-content-preview&box_client_version=2.14.1)

And here’s a picture of the inspiration (All credit to the show: “Is the order a rabbit?” for the image): 
![Inspiration](https://dl.boxcloud.com/api/2.0/internal_files/507675529988/versions/537556945988/representations/png_paged_2048x2048/content/1.png?access_token=1!MRrYQ-du844HnDsHT3t4FZI7nazKiBPFOJCEEAAcQIocL3KxFvCejKdfSOepps0VE72uCZnH3fPFU81hHdJR55D5murbzrLRQhh0-lXOGoX5c8bBMOzyQWJ_IGPrQp12munPy36czksnzZaYRjYXmjCkEr_EhrMkEk5KBlY-yayutnm7mHSFiKtPeDjqksrEzRMSKClYVUPw4nZoh3kpHNP8I2McBoPrBKmjbVyx9xQJcCxX3qVWR6tnOwD2Ayl7wPJw0udXyAWoDqooiVfw1fyDs7q8sBFDWdTU2hp6S4a5W7OScM0mhKbW9k_oeh-KJ5W07m6Licgkfyvt76T2mr-HBmPyQYqp-AxvBeXNgNb9YBvlun7OK6HeYOaxbQQZYk0bj1vKN_ugRcRh3VNpKmiSMlW9R2IC8U7VOgENg5jsj6RozDwwbXtJLp_t2IqF7Eb_BlT2Xap50myYKxpC2lwSNNNMzD3v2_em8qKC7ZsQ-6GwEKjyGCY1H5TiX2yVT6T4AGILXq-Wl4QBETfZdjqi6ivO3NsPedhG4rNcyF0fSX34VYFjNwDO8efi7i9tvw..&box_client_name=box-content-preview&box_client_version=2.14.1)

del *.class
del *.obj
cd ..
cd src\defaultPackage
javac *.java -Xlint:unchecked
cd ..
javah.exe -classpath . defaultPackage.JavaToC
copy defaultPackage_JavaToC.h "..\c&cuda"
cd ..
cd "c&cuda"
cl -c -Ox -I"C:\Program Files\Java\jdk1.7.0_25\include" -I"C:\Program Files\Java\jdk1.7.0_25\include\win32" Diffraction.c
nvcc -arch=sm_21 cuDiffraction.cu -c -m64 -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0\include" -I"C:\ProgramData\NVIDIA Corporation\CUDA Samples\v5.0\common\inc"
link /DLL /LIBPATH:"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0\lib\x64" /OUT:Diffraction64.dll *.obj cudart.lib cuda.lib /LIBPATH:"C:\ProgramData\NVIDIA Corporation\CUDA Samples\v5.0\common\lib"
del *.class
del *.obj
del "..\Diffraction64.dll"
copy Diffraction64.dll "..\"
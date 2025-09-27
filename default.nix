with import <nixpkgs> {};

let
in
mkShell {
    name = "cuda-and-openmp";
    stdenv = gcc12Stdenv;
    packages = [
        glm
        glfw
        freetype
        valgrind
        kdePackages.kcachegrind
    ];

    buildInputs = with pkgs; [
        (pkgs.stdenv.mkDerivation {
            pname = "sfml";
            version = "3.0.1";

            src = pkgs.fetchgit {
                url = "https://github.com/SFML/SFML.git";
                rev = "3.0.1";  # or a specific commit hash
                sha256 = "sha256-YqlrY0iIsxcjlLb+buMU0zpXo7/eKSKxOsITWf7BX6s=";
            };
            nativeBuildInputs = [ 
                cmake 
                pkg-config 
                xorg.libX11
                xorg.libXrandr
                xorg.libXinerama
                xorg.libXcursor
                xorg.libXi
                libGLU 
                libGL
                udev
                freetype
                libvorbis
                flac
            ];
        })
        pkgs.pkg-config 
        cudaPackages.cuda_cudart
        cudaPackages.cuda_nvcc
        cudaPackages.cuda_cccl
        cudaPackages.cudatoolkit
        cudaPackages.cuda_nvprof
        cudaPackages.cuda_nvvp
        jdk8
        linuxPackages.nvidia_x11
        gcc12
        libGLU libGL
        glm
        glfw
        freetype
        vulkan-loader
        pkg-config
        xorg.libX11
        xorg.libXrandr
        xorg.libXinerama
        xorg.libXcursor
        xorg.libXi
        cmake
        libvorbis
        flac
    ];

    shellHook = ''
        export SFML_PATH=${sfml}/lib/cmake
        export JAVA_HOME=${pkgs.jdk8}
        export CC=${pkgs.gcc12}/bin/gcc
        export CXX=${pkgs.gcc12}/bin/g++

        export CUDAHOSTCXX=${pkgs.gcc12}/bin/g++
        export CUDA_HOST_COMPILER=${pkgs.gcc12}/bin/gcc

        export CUDA_HOME=${pkgs.cudaPackages.cuda_cudart}
        export CUDA_PATH=${pkgs.cudaPackages.cuda_cudart}
        #
        export LD_LIBRARY_PATH=${pkgs.cudaPackages.cuda_cudart}/lib64:${pkgs.cudaPackages.cuda_cudart}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=${stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=${pkgs.linuxPackages.nvidia_x11}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH="$JAVA_HOME/lib:$JAVA_HOME/lib/server:$LD_LIBRARY_PATH"
        export LD_LIBRARY_PATH=${libvorbis}/lib:$LD_LIBRARY_PATH
        export LD_LIBRARY_PATH=${flac}/lib:$LD_LIBRARY_PATH

        export LIBRARY_PATH=${pkgs.cudaPackages.cuda_cudart}/lib64:${pkgs.cudaPackages.cuda_cudart}/lib:$LIBRARY_PATH
        export PATH=$SFML_PATH/bin:$PATH
    '';
}

{ pkgs ? import <nixpkgs> { } }:

let
  # nixpkgs only ships SFML 2.6, so we build 3.x from source.
  sfml = pkgs.stdenv.mkDerivation (finalAttrs: {
    pname = "sfml";
    version = "3.0.1";

    src = pkgs.fetchgit {
      url = "https://github.com/SFML/SFML.git";
      rev = "3.0.1"; # verify: `nix-prefetch-git https://github.com/SFML/SFML.git --rev 3.0.1`
      hash = "sha256-YqlrY0iIsxcjlLb+buMU0zpXo7/eKSKxOsITWf7BX6s=";
    };

    nativeBuildInputs = with pkgs; [
      cmake
      ninja
      pkg-config
    ];

    buildInputs = with pkgs; [
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
      openal          # SFML 3 audio backend
      vulkan-loader
    ];

    cmakeFlags = [
      "-DBUILD_SHARED_LIBS=ON"
      "-DSFML_BUILD_NETWORK=OFF"
      "-DSFML_BUILD_EXAMPLES=OFF"
      "-DSFML_BUILD_DOC=OFF"
    ];
  });
in
pkgs.mkShell {
  name = "cuda-and-openmp";

  packages = with pkgs; [
    glm
    glfw
    freetype
    valgrind
    kdePackages.kcachegrind
  ];

  buildInputs = with pkgs; [
    sfml
    pkg-config
    cmake

    # CUDA
    cudaPackages.cuda_cudart
    cudaPackages.cuda_nvcc
    cudaPackages.cuda_cccl
    cudaPackages.cudatoolkit
    cudaPackages.cuda_nvtx
    cudaPackages.cuda_gdb
    linuxPackages.nvidia_x11

    # JVM
    jdk8

    # Graphics / windowing
    libGLU
    libGL
    glm
    glfw
    vulkan-loader
    xorg.libX11
    xorg.libXrandr
    xorg.libXinerama
    xorg.libXcursor
    xorg.libXi

    # Audio / misc
    libvorbis
    flac
    openal
    freetype
  ];

  shellHook = ''
    export SFML_DIR=${sfml}/lib/cmake/SFML
    export JAVA_HOME=${pkgs.jdk8}

    export CUDAHOSTCXX=${pkgs.gcc}/bin/g++
    export CUDA_HOST_COMPILER=${pkgs.gcc}/bin/gcc
    export CUDA_HOME=${pkgs.cudaPackages.cuda_cudart}
    export CUDA_PATH=${pkgs.cudaPackages.cuda_cudart}

    export LD_LIBRARY_PATH=${pkgs.lib.makeLibraryPath [
      pkgs.cudaPackages.cuda_cudart
      pkgs.stdenv.cc.cc.lib
      pkgs.linuxPackages.nvidia_x11
      pkgs.libvorbis
      pkgs.flac
      pkgs.openal
      pkgs.libGL
      pkgs.glfw
      pkgs.vulkan-loader
      sfml
    ]}:$LD_LIBRARY_PATH

    export LIBRARY_PATH=${pkgs.lib.makeLibraryPath [
      pkgs.cudaPackages.cuda_cudart
    ]}:$LIBRARY_PATH

    export PATH=${sfml}/bin:$PATH
  '';
}

{
  description = "vFind development and build environment";

  inputs = {
    nixpkgs.url = "github:cachix/devenv-nixpkgs/rolling";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
    devenv.inputs.nixpkgs.follows = "nixpkgs";
    nixpkgs-python.url = "github:cachix/nixpkgs-python";
    nixpkgs-python.inputs.nixpkgs.follows = "nixpkgs";
    fenix.url = "github:nix-community/fenix";
    fenix.inputs.nixpkgs.follows = "nixpkgs";
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, devenv, systems, ... } @ inputs:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      packages = forEachSystem (system: {
        devenv-up = self.devShells.${system}.default.config.procfileScript;
      });

      devShells = forEachSystem
        (system:
          let
            pkgs = nixpkgs.legacyPackages.${system};
          in
          {
            default = devenv.lib.mkShell {
              inherit inputs pkgs;
              modules = [
                {
                  packages = with pkgs; [
                    cmake
                    libclang
                    pkg-config
                    (maturin.overrideAttrs (oldAttrs: rec {
                      version = "1.5.1",
                      src = pkgs.fetchFromGitHub {
                        owner = "PyO3";
                        repo = "maturin";
                        rev = "7d711f0c4a7c052608dc2e16d5c6721b9666d076";
                        hash = "sha256-3rID2epV1pCwpofFf9Wuafs1SlBWH7e7/4HPaSUAriQ=";
                        fetchSubmodules = true;
                      }
                    }))
                    # python3 packages
                    pkgs.python312Packages.pytest
                    pkgs.python312Packages.polars
                  ];

                  languages.python = {
                    enable = true;
                    version = "3.12";
                    manylinux.enable = system == "x86_64-linux" || system == "aarch64-linux";
                  };

                  languages.rust.enable = true;
                  languages.rust.channel = "nightly";
                  enterShell = ''
                    export LIBCLANG_PATH="$DEVENV_PROFILE/lib"
                  '';

                }
              ];
            };
          });
    };
}

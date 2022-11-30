{
  description = "ROOT dev env.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlay = (final: prev: {
          root-dev = prev.callPackage ./nix/default.nix {
            python = final.python3;
            inherit (prev.darwin.apple_sdk.frameworks) Cocoa CoreSymbolication OpenGL;
            noSplash = true;
          };
        });
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ overlay ];
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      {
        devShell = pkgs.mkShell {
          name = "root-dev-shell";
          buildInputs = with pythonPackages; [
            pkgs.clang-tools # For clang-format
            pkgs.root-dev

            # Linters
            pylint
          ];
        };
      });
}

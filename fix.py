import pathlib


def main():
    paths = list(pathlib.Path(".").glob("**/*.rg.rg"))

    print(f"Found {len(paths)} files to rename")
    print("Rename all files? (y/n)")
    answer = input()

    if answer != "y":
        print("Aborting")
        return

    print("Renaming files...")

    for path in paths:
        new_name = path.name.replace(".rg.rg", ".rg")

        if path.with_name(new_name).exists():
            print(f"File {path} already exists, skipping")
            continue

        path.rename(new_name)


if __name__ == "__main__":
    main()

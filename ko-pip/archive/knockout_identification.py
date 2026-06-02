def identify_knockouts(depth_file, knockout_file,non_knockout_file, percentage):
    # First pass to calculate average depth
    total_depth = 0
    num_positions = 0
    with open(depth_file, 'r') as depth:
        for line in depth:
            _, position, depth = line.strip().split('\t')
            position, depth = map(int, [position, depth])
            total_depth += depth
            num_positions += 1
    average_depth = total_depth / num_positions

    # Calculate threshold as percentage of average depth
    threshold = percentage / 100 * average_depth

    # Second pass to identify knockout regions
    knockout_regions = []
    non_deleted_regions = []
    start_ko = None
    end_last_ko = 0
    last_position = None
    with open(depth_file, 'r') as depth:
        for line in depth:
            _, position, depth = line.strip().split('\t')
            position, depth = map(int, [position, depth])
            if last_position is None:  # This is the first position
                last_position = position
            if depth < threshold:
                if start_ko is None:
                    start_ko = position
                if end_last_ko is not None and (start_ko - end_last_ko > 1):
                    non_deleted_regions.append((end_last_ko, start_ko - 1))
            elif start_ko is not None:
                knockout_regions.append((start_ko, position - 1))
                end_last_ko = position
                start_ko = None
            last_position = position

        if start_ko is not None:
            knockout_regions.append((start_ko, last_position))
        elif end_last_ko is not None and end_last_ko < last_position:
            non_deleted_regions.append((end_last_ko, last_position))
        elif end_last_ko is None:  # No knockout regions found
            non_deleted_regions.append((0, last_position))

    # Write output files
    with open(knockout_file, 'w') as out:
        for start, end in knockout_regions:
            out.write(f"{start}-{end}\n")

    with open(non_knockout_file, 'w') as out:
        # Remove duplicate entries
        non_deleted_regions = list(set(non_deleted_regions))
        # Sort regions
        non_deleted_regions.sort(key=lambda x: x[0])
        for start, end in non_deleted_regions:
            out.write(f"{start}\t{end}\n")

    print(f"[Takeda Pipeline Message] Knockout regions written to {knockout_file}")
    print(f"[Takeda Pipeline Message] Non-deleted regions written to {non_knockout_file}")


if __name__ == "__main__":
    pass
